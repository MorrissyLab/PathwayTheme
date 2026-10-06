"""Tests for the metadata-by-component test and the diagnostic plots."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pathwaytheme.pca.metadata import (associate_metadata, associate_levels,
                                       attribution_matrix, level_matrix)
from pathwaytheme.viz import (draw_workflow, draw_scree, draw_pca_scatter,
                              draw_theme_heatmap,
                              draw_metadata_attribution, draw_volcano,
                              draw_top_bars, draw_effect_concordance,
                              draw_level_attribution, render_scree,
                              render_theme_heatmap, render_level_attribution)


@pytest.fixture
def scores() -> pd.DataFrame:
    """PC1 separates group A from B by construction; PC2 is noise."""
    rng = np.random.default_rng(0)
    n = 30
    pc1 = np.concatenate([rng.normal(-3, 0.5, n), rng.normal(3, 0.5, n)])
    pc2 = rng.normal(0, 1, 2 * n)
    pc3 = rng.normal(0, 1, 2 * n)
    idx = [f"s{i}" for i in range(2 * n)]
    return pd.DataFrame({"PC1": pc1, "PC2": pc2, "PC3": pc3}, index=idx)


@pytest.fixture
def metadata(scores) -> pd.DataFrame:
    n = len(scores) // 2
    return pd.DataFrame({
        "group": ["A"] * n + ["B"] * n,
        "batch": (["x", "y"] * n)[: 2 * n],
        "score": np.linspace(0, 1, 2 * n),
        "constant": ["same"] * 2 * n,
        "mostly_missing": [None] * (2 * n - 2) + ["a", "b"],
    }, index=scores.index)


def test_finds_the_planted_association(scores, metadata):
    a = associate_metadata(scores, metadata)
    top = a.nlargest(1, "effect").iloc[0]
    assert top["component"] == "PC1"
    assert top["variable"] == "group"
    assert top["significant"]
    # PC2 is noise: it must not be attributed to the group, and the planted
    # effect must dwarf it.  A relative check, because the absolute eta^2 of a
    # rank test on perfectly separated groups is not 1 -- here it is ~0.75.
    pc2 = a[(a.component == "PC2") & (a.variable == "group")].iloc[0]
    assert not pc2["significant"]
    assert top["effect"] > 20 * max(pc2["effect"], 1e-3)


def test_reports_direction_and_test_used(scores, metadata):
    a = associate_metadata(scores, metadata)
    row = a[(a.component == "PC1") & (a.variable == "group")].iloc[0]
    assert row["test"] == "kruskal" and row["kind"] == "categorical"
    assert row["top_level"] == "B"       # B is the positive end of PC1
    cont = a[a.variable == "score"].iloc[0]
    assert cont["test"] == "spearman" and cont["effect_name"] == "rho_squared"


def test_skips_untestable_variables(scores, metadata):
    a = associate_metadata(scores, metadata, min_group_size=3)
    tested = set(a["variable"])
    assert "constant" not in tested          # one level
    assert "mostly_missing" not in tested    # every level below the floor
    assert {"group", "batch", "score"} <= tested


def test_fdr_is_over_the_whole_grid(scores, metadata):
    a = associate_metadata(scores, metadata)
    assert (a["fdr"] >= a["p_value"] - 1e-12).all()
    assert a["fdr"].le(1.0).all()
    # one row per (component, variable) pair, no duplicates
    assert not a.duplicated(["component", "variable"]).any()


def test_attribution_matrix_orders_components_numerically(scores, metadata):
    a = associate_metadata(scores, metadata)
    m = attribution_matrix(a)
    assert list(m.index) == ["PC1", "PC2", "PC3"]
    assert (m.fillna(0) >= 0).all().all()   # effects are non-negative
    signed = attribution_matrix(a, signed=True)
    assert signed.loc["PC1", "score"] != 0


def test_empty_inputs_do_not_raise(scores):
    assert associate_metadata(scores, pd.DataFrame(index=scores.index)).empty
    assert attribution_matrix(pd.DataFrame()).empty


# ── one-versus-rest: which level, not just which variable ──────────────────
@pytest.fixture
def one_level_driven():
    """PC1 separates C only; A and B are interchangeable on it.

    The omnibus test sees "subtype matters on PC1" without being able to say
    that only C is responsible -- which is what associate_levels reports.
    """
    rng = np.random.default_rng(7)
    n = 20
    pc1 = np.concatenate([rng.normal(0, 1, n), rng.normal(0, 1, n),
                          rng.normal(8, 1, n)])
    pc2 = rng.normal(0, 1, 3 * n)
    idx = [f"s{i}" for i in range(3 * n)]
    s = pd.DataFrame({"PC1": pc1, "PC2": pc2}, index=idx)
    m = pd.DataFrame({"subtype": ["A"] * n + ["B"] * n + ["C"] * n,
                      "sex": (["f", "m"] * (3 * n))[: 3 * n]}, index=idx)
    return s, m


def test_level_test_names_the_driving_level(one_level_driven):
    s, m = one_level_driven
    lv = associate_levels(s, m)
    pc1 = lv[(lv.component == "PC1") & (lv.variable == "subtype")]
    pc1 = pc1.set_index("level")
    # C is fully separated from the rest, and it is the level that ranks first
    assert pc1.loc["C", "significant"] and pc1.loc["C", "effect"] > 0.95
    assert pc1.loc["C", "auc"] > 0.97
    # A and B are also "significant", because the pooled rest they are compared
    # against contains C -- the ranking, not the flag, is what isolates C
    assert pc1["effect"].abs().idxmax() == "C"
    assert pc1.loc["C", "effect"] > 1.8 * pc1.loc[["A", "B"], "effect"].abs().max()
    # A and B are indistinguishable from each other
    assert pc1.loc["A", "effect"] == pytest.approx(pc1.loc["B", "effect"], abs=0.25)
    # and the omnibus test on the same data cannot make that distinction
    omni = associate_metadata(s, m)
    row = omni[(omni.component == "PC1") & (omni.variable == "subtype")].iloc[0]
    assert row["significant"] and row["effect"] > 0.5


def test_level_effect_is_signed(one_level_driven):
    s, m = one_level_driven
    lv = associate_levels(s, m)
    pc1 = lv[(lv.component == "PC1") & (lv.variable == "subtype")]
    pc1 = pc1.set_index("level")
    assert pc1.loc["C", "effect"] > 0                 # C sits high on PC1
    assert pc1.loc["C", "direction"] == "high"
    # A and B are dragged below the mean by C's presence in "rest"
    assert pc1.loc["A", "effect"] < 0
    assert pc1.loc["A", "mean_level"] < pc1.loc["A", "mean_rest"]
    # rank-biserial and AUC are two views of the same number
    assert np.allclose(pc1["effect"], 2 * pc1["auc"] - 1)


def test_binary_variable_is_not_double_counted(one_level_driven):
    s, m = one_level_driven
    lv = associate_levels(s, m)
    sex = lv[(lv.component == "PC1") & (lv.variable == "sex")]
    assert len(sex) == 2 and int(sex["redundant"].sum()) == 1
    # mirror rows: same p, same adjusted p, opposite sign
    assert sex["p_value"].nunique() == 1
    assert sex["fdr"].nunique() == 1
    assert sex["effect"].sum() == pytest.approx(0.0)
    assert lv["fdr"].notna().all()


def test_level_matrix_labels_and_order(one_level_driven):
    s, m = one_level_driven
    lv = associate_levels(s, m)
    allvars = level_matrix(lv)
    assert list(allvars.index) == ["PC1", "PC2"]
    assert "subtype: C" in allvars.columns and "sex: f" in allvars.columns
    one = level_matrix(lv, variable="subtype")
    assert set(one.columns) == {"A", "B", "C"}
    assert one.loc["PC1", "C"] > 0
    assert (level_matrix(lv, signed=False) >= 0).all().all()
    assert level_matrix(pd.DataFrame()).empty
    assert level_matrix(lv, variable="not_a_column").empty


def test_level_test_skips_continuous_and_tiny_levels(scores, metadata):
    lv = associate_levels(scores, metadata, min_group_size=3)
    tested = set(lv["variable"])
    assert "score" not in tested          # continuous
    assert "constant" not in tested       # one level
    assert "mostly_missing" not in tested
    assert {"group", "batch"} <= tested
    assert associate_levels(scores, pd.DataFrame(index=scores.index)).empty


# ── plots: they must draw onto a supplied ax, and survive empty input ──────
def _diff_table() -> pd.DataFrame:
    rng = np.random.default_rng(1)
    n = 200
    return pd.DataFrame({
        "comparison": ["A_vs_B"] * n,
        "pathway": [f"GOBP_TERM_{i}" for i in range(n)],
        "effect": rng.normal(0, 0.2, n),
        "fdr": rng.uniform(0, 1, n),
    })


def _category_summary() -> pd.DataFrame:
    return pd.DataFrame({
        "comparison": ["A_vs_B"] * 4 + ["C_vs_D"] * 4,
        "significance": ["significant"] * 8,
        "category": ["immune", "metabolism", "cell cycle", "unmapped"] * 2,
        "mean_effect": [0.2, -0.1, 0.05, 0.0, -0.2, 0.1, -0.05, 0.0],
        "n_pathways": [10, 5, 8, 99] * 2,
        "n_up": [10, 0, 6, 50] * 2,
        "n_down": [0, 5, 2, 49] * 2,
    })


def test_draw_functions_accept_an_ax(scores, metadata):
    var = np.array([0.5, 0.3, 0.2])
    a = associate_metadata(scores, metadata)
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    fl = axes.ravel()
    draw_scree(fl[0], var, elbow=2)
    draw_pca_scatter(fl[1], scores, var, colour_by=metadata["group"],
                     marker_by=metadata["batch"], group_means=True)
    draw_metadata_attribution(fl[2], a, variance_explained=var,
                              technical=["batch"])
    draw_theme_heatmap(fl[3], _category_summary())
    draw_volcano(fl[4], _diff_table(), comparison="A_vs_B")
    draw_top_bars(fl[5], _diff_table(), per_side=4)
    rho = draw_effect_concordance(fl[6], pd.Series([1.0, 2, 3]),
                                  pd.Series([1.0, 2, 3]))
    assert rho == pytest.approx(1.0)
    # every axis that was asked to draw actually has artists on it
    for ax in fl[:7]:
        assert ax.get_children()
    plt.close(fig)


def test_effect_concordance_says_which_way_and_names_the_points():
    """``poles`` writes the direction onto both axes; ``annotate`` names points.

    A scatter of two effect columns is unreadable without both: nothing on the
    axes says which end is which, and an anonymous cloud gives a reader no way
    to check the correlation against anything they know.  The disagreements are
    the points an agreement claim rests on, so they are named whatever their
    size.
    """
    fig, ax = plt.subplots()
    x = pd.Series([0.5, 0.3, -0.4, -0.2, 0.1], index=list("abcde"))
    y = pd.Series([0.4, 0.2, -0.5, -0.3, -0.1], index=list("abcde"))
    rho = draw_effect_concordance(
        ax, x, y, unit="themes", colour_by_agreement=True, annotate=2,
        annotate_disagreements=True, poles=("higher in B", "higher in A"))
    assert rho == pytest.approx(1.0)          # ranks agree; one sign does not
    # both axes carry the same quantity, so the direction is written once, in
    # the two corners where the analyses agree -- not twice onto the axes
    assert "higher in A" not in ax.get_xlabel()
    assert "higher in A" not in ax.get_ylabel()
    corners = {t.get_text() for t in ax.texts}
    assert {"higher in A", "higher in B"} <= corners
    # the subtitle is read by someone who has not read the caption, so it says
    # what the exceptions do, in words.  The count of shared points is the same
    # in every panel of a row, so it lives in the caption instead.
    assert "1 theme points the opposite way" in ax.get_title()
    assert "in common" not in ax.get_title()
    named = corners - {"higher in A", "higher in B"}
    assert named == {"a", "c", "e"}           # one each end, plus the flip
    # the identity-line label belongs to the unlabelled form, where nothing
    # else on the panel explains the diagonal
    plt.close(fig)

    fig, ax = plt.subplots()
    draw_effect_concordance(ax, x, y)
    assert "equal effect" in " ".join(t.get_text() for t in ax.texts)
    plt.close(fig)


def test_theme_heatmap_drops_unmapped_by_default():
    c = _category_summary()
    fig, ax = plt.subplots()
    draw_theme_heatmap(ax, c)
    assert "unmapped" not in [t.get_text() for t in ax.get_yticklabels()]
    plt.close(fig)
    fig, ax = plt.subplots()
    draw_theme_heatmap(ax, c, drop_unmapped=False)
    assert "unmapped" in [t.get_text() for t in ax.get_yticklabels()]
    plt.close(fig)


def test_plots_tolerate_empty_frames():
    empty_diff = pd.DataFrame(columns=["comparison", "pathway", "effect", "fdr"])
    fig, axes = plt.subplots(1, 4)
    draw_volcano(axes[0], empty_diff)
    draw_top_bars(axes[1], empty_diff)
    draw_theme_heatmap(axes[2], pd.DataFrame(
        columns=["comparison", "significance", "category", "mean_effect",
                 "n_pathways"]))
    draw_metadata_attribution(axes[3], pd.DataFrame())
    plt.close(fig)
    fig, ax = plt.subplots()
    draw_level_attribution(ax, pd.DataFrame())
    plt.close(fig)


def test_level_attribution_draws_and_renders(one_level_driven, tmp_path):
    s, m = one_level_driven
    lv = associate_levels(s, m)
    fig, ax = plt.subplots()
    draw_level_attribution(ax, lv, variable="subtype",
                           variance_explained=np.array([0.7, 0.3]))
    assert [t.get_text() for t in ax.get_yticklabels()] == ["PC1 (70%)", "PC2 (30%)"]
    assert ax.get_children()
    plt.close(fig)
    p = render_level_attribution(lv, tmp_path / "levels.pdf")
    assert p.exists() and p.stat().st_size > 0


def test_render_wrappers_write_files(tmp_path):
    p = render_scree(np.array([0.6, 0.3, 0.1]), tmp_path / "scree.pdf", elbow=2)
    assert p.exists() and p.stat().st_size > 0
    p = render_theme_heatmap(_category_summary(), tmp_path / "themes.pdf")
    assert p.exists() and p.stat().st_size > 0


# ── real metadata has gaps, and real cohorts are large ────────────────────
def test_scatter_skips_the_mean_of_a_tiny_group(one_level_driven):
    """A group below the floor gets points but no mean line.

    A dashed line at a single sample reads as a group position; the cohort this
    was found on has an 8-sample and a 1-sample level next to a 658-sample one.
    """
    s, m = one_level_driven
    tiny = m["subtype"].copy()
    tiny.iloc[0] = "solo"                     # one sample in its own level
    fig, ax = plt.subplots()
    draw_pca_scatter(ax, s, np.array([0.7, 0.3]), colour_by=tiny,
                     group_means=True, min_group_for_mean=3)
    means = [ln for ln in ax.get_lines() if ln.get_linestyle() == "--"]
    assert len(means) == 3                    # A, B, C -- not "solo"
    plt.close(fig)


def test_scatter_can_draw_any_component_pair(one_level_driven):
    s, m = one_level_driven
    fig, ax = plt.subplots()
    draw_pca_scatter(ax, s, np.array([0.7, 0.3]), components=("PC2", "PC1"),
                     colour_by=m["subtype"])
    assert ax.get_xlabel().startswith("PC2") and ax.get_ylabel().startswith("PC1")
    plt.close(fig)


def test_attribution_plots_take_a_horizontal_colourbar(one_level_driven):
    s, m = one_level_driven
    a = associate_metadata(s, m)
    lv = associate_levels(s, m)
    fig, axes = plt.subplots(1, 2)
    draw_metadata_attribution(axes[0], a, colorbar_horizontal=True)
    draw_level_attribution(axes[1], lv, colorbar_horizontal=True)
    # a horizontal bar is wider than it is tall; a vertical one is not
    bars = [ax for ax in fig.axes if ax not in list(axes)]
    assert len(bars) == 2
    for cb in bars:
        w, h = cb.get_position().width, cb.get_position().height
        assert w > h
    plt.close(fig)


def test_theme_heatmap_wraps_long_names():
    c = _category_summary()
    c["category"] = c["category"].replace(
        {"metabolism": "a very long theme name that would reach into the panel"})
    fig, ax = plt.subplots()
    draw_theme_heatmap(ax, c, label_chars=12)
    wrapped = [t.get_text() for t in ax.get_yticklabels() if "\n" in t.get_text()]
    assert wrapped, "a long theme name should be wrapped, not left to overflow"
    plt.close(fig)


@pytest.mark.parametrize("orientation, size",
                         [("horizontal", (7.3, 1.5)), ("vertical", (3.0, 4.6))])
def test_workflow_text_never_collides(orientation, size):
    """No label in the schematic may be drawn over another.

    The boxes used to be laid out in axes fractions with hand-set constants, so
    the drawing was only correct at the panel shape it had been tuned on: a
    method line that wrapped to two ran out through the bottom border, and the
    number badge sat on top of the first one.  Both failures show up here as
    two pieces of text sharing pixels, which is cheap to assert and is the
    thing a reader actually notices.
    """
    fig = plt.figure(figsize=size)
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
    draw_workflow(ax, orientation=orientation)
    r = fig.canvas.get_renderer()
    boxes = [(t.get_text(), t.get_window_extent(r)) for t in ax.texts]
    assert len(boxes) > 12, "the schematic should carry every stage's text"
    for i, (a_text, a) in enumerate(boxes):
        for b_text, b in boxes[i + 1:]:
            assert not a.overlaps(b), f"{a_text!r} is drawn over {b_text!r}"
    plt.close(fig)


def test_workflow_text_stays_inside_its_box():
    """Text that starts inside a stage box has to finish inside it."""
    from matplotlib.patches import FancyBboxPatch

    fig = plt.figure(figsize=(7.3, 1.5))
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
    draw_workflow(ax)
    r = fig.canvas.get_renderer()
    stages = [p.get_window_extent(r) for p in ax.patches
              if isinstance(p, FancyBboxPatch)]
    assert len(stages) == 6
    inside = 0
    for t in ax.texts:
        e = t.get_window_extent(r)
        for s in stages:
            if s.contains(*e.get_points()[0]):
                inside += 1
                assert s.contains(*e.get_points()[1]), \
                    f"{t.get_text()!r} starts inside a box and ends outside it"
                break
    assert inside >= 12, "each stage box holds a name and a method line"
    plt.close(fig)

