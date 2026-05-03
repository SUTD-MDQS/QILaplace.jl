# zT Circuit Compression Animation
# Animates the z-Transform circuit (QFT ∘ DT) compression into an MPO.
# 8 wires: n1,n1',n2,n2',n3,n3',n4,n4' (main + copy registers, interleaved)
# Three phases: DT compression, Copy register, QFT compression.
from manim import *
import shutil
import os
from pathlib import Path
import numpy as np

config.background_color = "#222831"

class DTCircuitCompression(Scene):
    def construct(self):
        Text.set_default(color="#DDDDDD")
        MathTex.set_default(color="#DDDDDD")
        Tex.set_default(color="#DDDDDD")

        # ================================================================
        # THEME COLORS
        # ================================================================
        COLOR_PURPLE = "#8E47B7"       # V tensor
        COLOR_GREEN = "#49989B"        # T (orthogonal center)
        COLOR_LIGHT_BLUE = "#297DB5"   # Compressed bonds & S
        COLOR_ORANGE_RED = "#B14628"   # Heavy uncompressed bonds
        COLOR_BLACK = BLACK
        COLOR_T_US = "#407274"         # Sky Blue for T and U·S
        COLOR_H = "#3CB371"            # H gate fill (unitary, QFT section)
        COLOR_BROWN_1 = "#D98E43"
        COLOR_BROWN_2 = "#A47243"
        COLOR_BROWN_3 = "#6A4329"
        COLOR_YELLOW = "#EAD54B"       # U tensor
        COLOR_U = "#E7A84A"           # U gate hatch color (distinct from R gates)

        COLOR_HD = COLOR_H           # Hd gate fill matches H gate color
        COLOR_R_DT = "#C6A86C"         # R gate fill in DT section
        COLOR_R_COPY = "#FA98D3"       # R gate fill in copy section

        # ================================================================
        # GATE AND WIRE DIMENSIONS
        # ================================================================
        GATE_WIDTH = 0.50
        WIRE_THICKNESS = 3.5
        BOND_HEAVY = 12
        BOND_LIGHT = 6
        MIN_GATE_GAP = 0.12
        GATE_CENTER_SPACING = GATE_WIDTH + MIN_GATE_GAP
        GATE_LABEL_SCALE = 0.5
        WIRE_LABEL_SCALE = 0.56
        HATCH_LINES_SQUARE_GATE = 20

        # ================================================================
        # WIRE LAYOUT
        # ================================================================
        n_main = 4
        n_total = 2 * n_main  # 8 wires
        WIRE_SPACING_MAIN_COPY = 0.62
        WIRE_SPACING_PAIR = 0.88
        CIRCUIT_CENTER_X = 0.0

        wire_y_positions = []
        y_cursor = 3.0
        for i in range(n_main):
            wire_y_positions.append(y_cursor)             # n_i (main)
            y_cursor -= WIRE_SPACING_MAIN_COPY
            wire_y_positions.append(y_cursor)             # n_i' (copy)
            if i < n_main - 1:
                y_cursor -= WIRE_SPACING_PAIR

        # Vertically center the full interleaved wire stack around y=0.
        y_top = max(wire_y_positions)
        y_bottom = min(wire_y_positions)
        y_center = 0.5 * (y_top + y_bottom)
        wire_y_positions = [y - y_center for y in wire_y_positions]

        # ================================================================
        # VISUAL PRIMITIVE CLASSES
        # ================================================================

        class QuantumWire(VGroup):
            def __init__(self, y_pos, length, start_x, idx, is_copy=False):
                super().__init__()
                self.y_pos = y_pos
                self.idx = idx
                self.is_copy = is_copy

                dash_args = {"dash_length": 0.12, "dashed_ratio": 0.5} if is_copy else {}
                if is_copy:
                    self.line = DashedLine(
                        start=RIGHT * start_x + UP * y_pos,
                        end=RIGHT * (start_x + length) + UP * y_pos,
                        stroke_width=WIRE_THICKNESS * 0.8,
                        color="#AAAAAA",
                        **dash_args,
                    )
                else:
                    self.line = Line(
                        start=RIGHT * start_x + UP * y_pos,
                        end=RIGHT * (start_x + length) + UP * y_pos,
                        stroke_width=WIRE_THICKNESS,
                        color=WHITE,
                    )

                circ_rad = (WIRE_THICKNESS * 1.0) / 100
                ep_color = "#AAAAAA" if is_copy else WHITE
                self.start_circ = Circle(radius=circ_rad, color=ep_color, fill_opacity=1).move_to(self.line.get_start())
                self.end_circ = Circle(radius=circ_rad, color=ep_color, fill_opacity=1).move_to(self.line.get_end())

                prime = "'" if is_copy else ""
                main_idx = (idx // 2) + 1
                self.label = MathTex(f"|n_{{{main_idx}{prime}}}\\rangle").scale(WIRE_LABEL_SCALE).next_to(self.start_circ, LEFT, buff=0.10)
                self.add(self.line, self.start_circ, self.end_circ, self.label)
                self.set_z_index(-2) # Move wires behind everything

        class QuantumGate(VGroup):
            def __init__(self, label_str, x_pos, y_pos, color):
                super().__init__()
                self.x_pos = x_pos
                self.target_y = y_pos
                side = GATE_WIDTH
                self.box = RoundedRectangle(
                    corner_radius=0.06, width=side, height=side,
                    color=COLOR_BLACK, stroke_width=WIRE_THICKNESS * side,
                    fill_color=color, fill_opacity=1,
                ).move_to(RIGHT * x_pos + UP * y_pos)
                self.label = MathTex(label_str, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE).move_to(self.box.get_center())
                self.add(self.box, self.label)

        class HatchedGate(VGroup):
            def __init__(self, label_str, x_pos, y_pos, color, label_scale=GATE_LABEL_SCALE):
                super().__init__()
                self.x_pos = x_pos
                self.target_y = y_pos
                side = GATE_WIDTH
                self.bg = RoundedRectangle(
                    corner_radius=0.06, width=side, height=side,
                    color=COLOR_BLACK, stroke_width=0,
                    fill_color="#D3D3D3", fill_opacity=1.0,
                ).move_to(RIGHT * x_pos + UP * y_pos)
                hatch_group = VGroup()
                n_lines = HATCH_LINES_SQUARE_GATE
                step = (2 * side) / (n_lines + 1)
                for i in range(1, n_lines + 1):
                    c = -side + i * step
                    x1 = max(-side/2, c - side/2)
                    y1 = c - x1
                    x2 = min(side/2, c + side/2)
                    y2 = c - x2
                    line = Line(
                        RIGHT * (x_pos + x1) + UP * (y_pos + y1),
                        RIGHT * (x_pos + x2) + UP * (y_pos + y2),
                        stroke_width=2.0, color=color, stroke_opacity=1.0
                    )
                    hatch_group.add(line)
                self.box = RoundedRectangle(
                    corner_radius=0.06, width=side, height=side,
                    color=COLOR_BLACK, stroke_width=WIRE_THICKNESS * side,
                    fill_opacity=0,
                ).move_to(RIGHT * x_pos + UP * y_pos)
                self.hatching = hatch_group
                self.label = MathTex(label_str, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(label_scale).move_to(self.box.get_center())

                self.bg.set_z_index(1)
                self.hatching.set_z_index(2)
                self.box.set_z_index(3)
                self.label.set_z_index(4)

                self.add(self.bg, self.hatching, self.box, self.label)

        class InvertedPhaseTrain(VGroup):
            def __init__(self, x_pos, control_y, targets_info):
                super().__init__()
                self.targets_info = targets_info
                self.dot = Dot(point=RIGHT * x_pos + UP * control_y, radius=(WIRE_THICKNESS * 2.0)/100, color=WHITE)
                max_y = max([t[0] for t in targets_info])
                self.cline = Line(self.dot.get_center(), RIGHT * x_pos + UP * max_y, stroke_width=WIRE_THICKNESS, color=WHITE).set_z_index(-1)
                self.add(self.dot, self.cline)
                self.gates = []
                for (y_pos, label_str, color, is_hatched) in targets_info:
                    gate = HatchedGate(label_str, x_pos, y_pos, color) if is_hatched else QuantumGate(label_str, x_pos, y_pos, color)
                    gate.set_z_index(1, family=False)
                    self.gates.append(gate)
                    self.add(gate)

        class PhaseTrain(VGroup):
            def __init__(self, x_pos, control_y, targets_info):
                super().__init__()
                self.targets_info = targets_info
                self.dot = Dot(point=RIGHT * x_pos + UP * control_y, radius=(WIRE_THICKNESS * 2.0)/100, color=WHITE)
                min_y = min([t[0] for t in targets_info])
                self.cline = Line(self.dot.get_center(), RIGHT * x_pos + UP * min_y, stroke_width=WIRE_THICKNESS, color=WHITE).set_z_index(-1)
                self.add(self.dot, self.cline)
                self.gates = []
                for (y_pos, label_str, color, is_hatched) in targets_info:
                    gate = HatchedGate(label_str, x_pos, y_pos, color) if is_hatched else QuantumGate(label_str, x_pos, y_pos, color)
                    gate.set_z_index(1, family=False)
                    self.gates.append(gate)
                    self.add(gate)

        class TensorNode(VGroup):
            def __init__(self, label_str, x_pos, y_pos, t_type="U"):
                super().__init__()
                color = COLOR_YELLOW if t_type == "U" else (COLOR_PURPLE if t_type == "V" else COLOR_GREEN)
                self.box = RoundedRectangle(
                    corner_radius=0.06, width=GATE_WIDTH, height=GATE_WIDTH,
                    color=COLOR_BLACK, stroke_width=WIRE_THICKNESS * GATE_WIDTH,
                    fill_color=color, fill_opacity=1,
                ).move_to(RIGHT * x_pos + UP * y_pos)
                self.label = MathTex(label_str, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE).move_to(self.box.get_center())
                self.add(self.box, self.label)
                self.set_z_index(5)

        # Factory helpers
        def make_wide_tensor(x_left, x_right, y_pos, label_str, color):
            width = x_right - x_left
            center_x = (x_left + x_right) / 2
            box = RoundedRectangle(
                corner_radius=0.06, width=width, height=GATE_WIDTH,
                color=COLOR_BLACK, stroke_width=WIRE_THICKNESS,
                fill_color=color, fill_opacity=1,
            ).move_to(RIGHT * center_x + UP * y_pos)
            label = MathTex(label_str, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE).move_to(box.get_center())
            group = VGroup(box, label)
            group.set_z_index(5)
            return group

        def make_tall_tensor(x_pos, y_top, y_bottom, label_str, color, level_padding=0.3):
            height = abs(y_top - y_bottom) + (2 * level_padding)
            center_y = (y_top + y_bottom) / 2
            box = RoundedRectangle(
                corner_radius=0.06, width=GATE_WIDTH, height=height,
                color=COLOR_BLACK, stroke_width=WIRE_THICKNESS,
                fill_color=color, fill_opacity=1,
            ).move_to(RIGHT * x_pos + UP * center_y)
            label = MathTex(label_str, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE).move_to(box.get_center())
            group = VGroup(box, label)
            group.set_z_index(5)
            return group

        def make_wide_hatched_tensor(x_left, x_right, y_pos, color):
            width = x_right - x_left
            center_x = (x_left + x_right) / 2
            half_w = width / 2
            half_h = GATE_WIDTH / 2

            bg = RoundedRectangle(
                corner_radius=0.06, width=width, height=GATE_WIDTH,
                color=COLOR_BLACK, stroke_width=0,
                fill_color="#D3D3D3", fill_opacity=1.0,
            ).move_to(RIGHT * center_x + UP * y_pos)

            hatch_group = VGroup()
            c_min = -half_w - half_h
            c_max = half_w + half_h
            square_pitch = (2 * GATE_WIDTH) / (HATCH_LINES_SQUARE_GATE + 1)
            n_lines = max(1, int(round((c_max - c_min) / square_pitch)) - 1)
            step = (c_max - c_min) / (n_lines + 1)
            for i in range(1, n_lines + 1):
                c = c_min + i * step
                x1 = max(-half_w, c - half_h)
                y1 = c - x1
                x2 = min(half_w, c + half_h)
                y2 = c - x2
                line = Line(
                    RIGHT * (center_x + x1) + UP * (y_pos + y1),
                    RIGHT * (center_x + x2) + UP * (y_pos + y2),
                    stroke_width=1.6, color=color, stroke_opacity=1.0,
                )
                hatch_group.add(line)

            border = RoundedRectangle(
                corner_radius=0.06, width=width, height=GATE_WIDTH,
                color=COLOR_BLACK, stroke_width=WIRE_THICKNESS,
                fill_opacity=0,
            ).move_to(RIGHT * center_x + UP * y_pos)

            # Match z-index ordering used by HatchedGate (bg=1, hatching=2, border=3)
            bg.set_z_index(1)
            hatch_group.set_z_index(2)
            border.set_z_index(3)
            
            return VGroup(bg, hatch_group, border)

        def make_tall_hatched_tensor(x_pos, y_top, y_bottom, color):
            height = abs(y_top - y_bottom)
            center_y = (y_top + y_bottom) / 2
            half_w = GATE_WIDTH / 2
            half_h = height / 2

            bg = RoundedRectangle(
                corner_radius=0.06, width=GATE_WIDTH, height=height,
                color=COLOR_BLACK, stroke_width=0,
                fill_color="#D3D3D3", fill_opacity=1.0,
            ).move_to(RIGHT * x_pos + UP * center_y)

            hatch_group = VGroup()
            c_min = -half_w - half_h
            c_max = half_w + half_h
            square_pitch = (2 * GATE_WIDTH) / (HATCH_LINES_SQUARE_GATE + 1)
            n_lines = max(1, int(round((c_max - c_min) / square_pitch)) - 1)
            step = (c_max - c_min) / (n_lines + 1)
            for i in range(1, n_lines + 1):
                c = c_min + i * step
                x1 = max(-half_w, c - half_h)
                y1 = c - x1
                x2 = min(half_w, c + half_h)
                y2 = c - x2
                line = Line(
                    RIGHT * (x_pos + x1) + UP * (center_y + y1),
                    RIGHT * (x_pos + x2) + UP * (center_y + y2),
                    stroke_width=1.6, color=color, stroke_opacity=1.0,
                )
                hatch_group.add(line)

            border = RoundedRectangle(
                corner_radius=0.06, width=GATE_WIDTH, height=height,
                color=COLOR_BLACK, stroke_width=WIRE_THICKNESS,
                fill_opacity=0,
            ).move_to(RIGHT * x_pos + UP * center_y)

            # Match z-index ordering used by HatchedGate (bg=1, hatching=2, border=3)
            bg.set_z_index(1)
            hatch_group.set_z_index(2)
            border.set_z_index(3)
            
            return VGroup(bg, hatch_group, border)

        def make_square_hatched_tensor(x_pos, y_pos, color):
            return make_wide_hatched_tensor(
                x_pos - 0.5 * GATE_WIDTH,
                x_pos + 0.5 * GATE_WIDTH,
                y_pos,
                color,
            )

        def _boundary_anchor_for_point(box_mob, point):
            center = box_mob.get_center()
            cx, cy = center[0], center[1]
            dx, dy = point[0] - cx, point[1] - cy
            half_w = box_mob.width / 2
            half_h = box_mob.height / 2

            scale_x = abs(dx) / max(half_w, 1e-8)
            scale_y = abs(dy) / max(half_h, 1e-8)
            scale = max(scale_x, scale_y, 1e-8)
            return np.array([cx + dx / scale, cy + dy / scale, 0.0])

        def get_dt_wire_target_gates(dt_items, wire_y, train_indices, y_tol=1e-3):
            selected = []
            for train_idx in train_indices:
                if train_idx < 0 or train_idx >= len(dt_items):
                    continue
                train = dt_items[train_idx]
                if not isinstance(train, InvertedPhaseTrain):
                    continue
                for gate, target_info in zip(train.gates, train.targets_info):
                    if abs(target_info[0] - wire_y) < y_tol:
                        selected.append(gate)
                        break
            selected.sort(key=lambda gate: gate.get_center()[0])
            return selected

        def merge_wire_tt_block(wire_y, hd_gate, target_gates, ti_label, ti_color, bond_anchor_registry=None, show_label=True, remove_lines=None, attached_leg_lines=None):
            if not target_gates and hd_gate is None:
                raise ValueError("merge_wire_tt_block needs at least one source gate")
            remove_lines = remove_lines or []
            attached_leg_lines = attached_leg_lines or []

            consumed_gates = list(target_gates)
            removed_label_ids = set()

            def _x_bounds_for_source(src_gate):
                if isinstance(src_gate, Dot):
                    c = src_gate.get_center()[0]
                    return c - (GATE_WIDTH / 2), c + (GATE_WIDTH / 2)
                return src_gate.get_left()[0], src_gate.get_right()[0]

            if hd_gate is not None and target_gates:
                hd_x = hd_gate.get_center()[0]
                right_targets = [gate for gate in target_gates if gate.get_center()[0] >= hd_x]
                nearest = min(right_targets if right_targets else target_gates, key=lambda gate: abs(gate.get_center()[0] - hd_x))
                hd_gate.set_z_index(nearest.get_z_index() - 1)
                if hasattr(nearest, "label"):
                    self.play(
                        hd_gate.animate.move_to(nearest.get_center()),
                        Unwrite(nearest.label),
                        run_time=0.3,
                    )
                    removed_label_ids.add(id(nearest.label))
                else:
                    self.play(hd_gate.animate.move_to(nearest.get_center()), run_time=0.3)
                consumed_gates.insert(0, hd_gate)
            elif hd_gate is not None:
                consumed_gates.insert(0, hd_gate)

            x_bounds = [_x_bounds_for_source(gate) for gate in consumed_gates]
            x_left = min(lb for lb, _ in x_bounds) # no padding, we want the merged tensor to tightly cover the source gates
            x_right = max(rb for _, rb in x_bounds) # no padding, we want the merged tensor to tightly cover the source gates
            merged_shape = make_wide_hatched_tensor(x_left, x_right, wire_y, ti_color)
            ti_text = None
            if show_label:
                ti_text = MathTex(ti_label, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
                ti_text.move_to(merged_shape[2].get_center())
                ti_text.set_z_index(7)

            # Add updaters to attached legs so they track the gates as they merge.
            # We assume legs are vertical.
            for line in attached_leg_lines:
                if line.get_num_points() == 0: continue
                # Identify if the line is attached via top or bottom
                l_start, l_end = line.get_start(), line.get_end()
                is_top_attached = abs(l_end[1] - wire_y) < 0.1
                is_bottom_attached = abs(l_start[1] - wire_y) < 0.1
                
                if is_top_attached:
                    line.add_updater(lambda mob, l=line: mob.put_start_and_end_on(l.get_start(), [mob.get_end()[0], merged_shape[2].get_top()[1], 0]))
                elif is_bottom_attached:
                    line.add_updater(lambda mob, l=line: mob.put_start_and_end_on([mob.get_start()[0], merged_shape[2].get_bottom()[1], 0], l.get_end()))

            label_unwrites = [
                Unwrite(gate.label)
                for gate in consumed_gates
                if hasattr(gate, "label") and id(gate.label) not in removed_label_ids
            ]
            if label_unwrites:
                self.play(*label_unwrites, run_time=0.2)

            merge_targets = [merged_shape.copy() for _ in consumed_gates]
            merge_anims = []
            for gate, target_copy in zip(consumed_gates, merge_targets):
                if isinstance(gate, HatchedGate):
                    gate_shape = VGroup(gate.bg, gate.hatching, gate.box)
                elif isinstance(gate, QuantumGate):
                    gate_shape = gate.box
                else:
                    gate_shape = gate
                merge_anims.append(ReplacementTransform(gate_shape, target_copy))

            merge_anims.extend([FadeOut(line) for line in remove_lines])

            if show_label and ti_text is not None:
                self.play(*merge_anims, Write(ti_text), run_time=0.45)
            else:
                self.play(*merge_anims, run_time=0.45)
            
            for line in attached_leg_lines:
                line.clear_updaters()

            for target_copy in merge_targets:
                self.remove(target_copy)

            if show_label and ti_text is not None:
                self.add(merged_shape, ti_text)
            else:
                self.add(merged_shape)
            self.remove(*consumed_gates)

            anchor_points = {
                "left": merged_shape[2].get_left(),
                "right": merged_shape[2].get_right(),
                "top": merged_shape[2].get_top(),
                "bottom": merged_shape[2].get_bottom(),
            }

            source_anchor_map = {}
            for gate in consumed_gates:
                anchor = _boundary_anchor_for_point(merged_shape[2], gate.get_center())
                source_anchor_map[id(gate)] = anchor
                if bond_anchor_registry is not None:
                    bond_anchor_registry[id(gate)] = anchor

            return {
                "tensor": VGroup(merged_shape, ti_text) if (show_label and ti_text is not None) else VGroup(merged_shape),
                "shape": merged_shape,
                "label": ti_text,
                "anchors": anchor_points,
                "source_anchor_map": source_anchor_map,
                "consumed": consumed_gates,
                "wire_y": wire_y,
            }

        def split_tensor_to_right_canon(
            tensor_info,
            current_wire_y,
            neighbor_wire_y,
            right_canon_label,
            remainder_label,
            attached_leg_lines=None,
            right_canon_color=COLOR_PURPLE,
            remainder_color=COLOR_T_US,
            move_time=0.7,
        ):
            attached_leg_lines = attached_leg_lines or []

            ti_shape = tensor_info["shape"]
            ti_border = ti_shape[2]
            ti_label_mob = tensor_info["label"]
            ti_center = ti_border.get_center()
            remainder_center_x = ti_center[0]

            # Keep existing attached legs smoothly connected while remainder moves.
            for line in attached_leg_lines:
                start_point = line.get_start().copy()
                end_x = line.get_end()[0]

                def leg_updater(mob, src=start_point, x_fixed=end_x):
                    target = np.array([x_fixed, ti_shape[2].get_top()[1], 0.0])
                    mob.put_start_and_end_on(src, target)

                line.add_updater(leg_updater)

            # Step 1: Move tensor to midpoint while relabeling and creating right-canon tensor.
            midpoint_y = 0.5 * (current_wire_y + neighbor_wire_y)
            ti_target = np.array([remainder_center_x, midpoint_y, 0.0])
            remainder_label_mob = MathTex(remainder_label, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
            remainder_label_mob.move_to(ti_shape[2].get_center())
            remainder_label_mob.set_z_index(8)
            remainder_label_mob.add_updater(lambda mob: mob.move_to(ti_shape[2].get_center()))
            if ti_label_mob is not None:
                ti_label_mob.add_updater(lambda mob: mob.move_to(ti_shape[2].get_center()))

            right_canon_shape = make_wide_hatched_tensor(
                remainder_center_x - 0.5 * GATE_WIDTH,
                remainder_center_x + 0.5 * GATE_WIDTH,
                current_wire_y,
                right_canon_color,
            )
            right_canon_label_mob = MathTex(right_canon_label, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
            right_canon_label_mob.move_to(right_canon_shape[2].get_center())
            right_canon_label_mob.set_z_index(8)

            svd_bond = Line(
                right_canon_shape[2].get_center(),
                ti_shape[2].get_top(),
                stroke_width=BOND_HEAVY,
                color=COLOR_ORANGE_RED,
            ).set_z_index(-5)
            svd_bond.add_updater(
                lambda mob: mob.put_start_and_end_on(
                    right_canon_shape[2].get_center(),
                    ti_shape[2].get_top(),
                )
            )

            right_canon_sweep = AnimationGroup(
                GrowFromEdge(right_canon_shape[0], LEFT),
                LaggedStart(*[GrowFromEdge(h, LEFT) for h in right_canon_shape[1]], lag_ratio=0.03),
                GrowFromEdge(right_canon_shape[2], LEFT),
                Write(right_canon_label_mob),
                lag_ratio=0.0,
            )

            text_swap_anims = []
            if ti_label_mob is not None:
                text_swap_anims.append(Unwrite(ti_label_mob, run_time=move_time / 2))
            text_swap_anims.append(Write(remainder_label_mob, run_time=move_time / 2))
            text_swap = Succession(*text_swap_anims)

            self.play(
                ti_shape.animate.move_to(ti_target),
                text_swap,
                right_canon_sweep,
                Create(svd_bond),
                run_time=move_time,
            )

            if ti_label_mob is not None:
                self.remove(ti_label_mob)
                ti_label_mob.clear_updaters()
            remainder_label_mob.move_to(ti_shape[2].get_center())
            remainder_label_mob.clear_updaters()
            svd_bond.clear_updaters()
            svd_bond.put_start_and_end_on(right_canon_shape[2].get_center(), ti_shape[2].get_top())
            for line in attached_leg_lines:
                line.clear_updaters()
                start_point = line.get_start().copy()
                end_x = line.get_end()[0]
                line.put_start_and_end_on(start_point, np.array([end_x, ti_shape[2].get_top()[1], 0.0]))

            ti_shape[0].set_fill("#D3D3D3", opacity=1.0)
            for hatch in ti_shape[1]:
                hatch.set_color(remainder_color)

            self.add(right_canon_shape, right_canon_label_mob, svd_bond)

            result = {
                "right_canon_tensor": VGroup(right_canon_shape, right_canon_label_mob),
                "remainder_tensor": VGroup(ti_shape, remainder_label_mob),
                "svd_bond": svd_bond,
                "legs": VGroup(*attached_leg_lines),
                "right_canon_anchors": {
                    "left": right_canon_shape[2].get_left(),
                    "right": right_canon_shape[2].get_right(),
                    "top": right_canon_shape[2].get_top(),
                    "bottom": right_canon_shape[2].get_bottom(),
                    "center": right_canon_shape[2].get_center(),
                },
                "remainder_anchors": {
                    "left": ti_shape[2].get_left(),
                    "right": ti_shape[2].get_right(),
                    "top": ti_shape[2].get_top(),
                    "bottom": ti_shape[2].get_bottom(),
                    "center": ti_shape[2].get_center(),
                },
                "current_wire_y": current_wire_y,
                "neighbor_wire_y": neighbor_wire_y,
            }
            return result

        def split_tensor_to_left_canon(
            tensor_info,
            current_wire_y,
            neighbor_wire_y,
            left_canon_label,
            remainder_label,
            attached_leg_lines=None,
            left_canon_color=COLOR_U,
            remainder_color=COLOR_T_US,
            move_time=0.7,
        ):
            attached_leg_lines = attached_leg_lines or []

            ti_shape = tensor_info["shape"]
            ti_border = ti_shape[2]
            ti_label_mob = tensor_info["label"]
            ti_center = ti_border.get_center()
            remainder_center_x = ti_center[0]

            # Shift legs along with the remainder tensor vertically
            midpoint_y = 0.5 * (current_wire_y + neighbor_wire_y)
            y_shift = midpoint_y - current_wire_y
            
            for line in attached_leg_lines:
                start_point = line.get_start().copy()
                end_point = line.get_end().copy()
                def leg_updater(mob, s=start_point, e=end_point, sh=y_shift, t=ti_shape):
                    # Instead of stretching, we shift the entire line based on the tensor's motion relative to start
                    curr_sh = (t[2].get_center()[1] - current_wire_y)
                    mob.put_start_and_end_on(s + UP * curr_sh, e + UP * curr_sh)
                line.add_updater(leg_updater)

            ti_target = np.array([remainder_center_x, midpoint_y, 0.0])
            
            remainder_label_mob = MathTex(remainder_label, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
            remainder_label_mob.move_to(ti_shape[2].get_center())
            remainder_label_mob.set_z_index(8)
            remainder_label_mob.add_updater(lambda mob: mob.move_to(ti_shape[2].get_center()))
            if ti_label_mob is not None:
                ti_label_mob.add_updater(lambda mob: mob.move_to(ti_shape[2].get_center()))

            # New left-canon tensor at current_wire_y exactly centered at the site center
            left_canon_shape = make_square_hatched_tensor(
                remainder_center_x,
                current_wire_y,
                left_canon_color,
            )
            left_canon_label_mob = MathTex(left_canon_label, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
            left_canon_label_mob.move_to(left_canon_shape[2].get_center())
            left_canon_label_mob.set_z_index(8)

            svd_bond = Line(
                ti_shape[2].get_bottom(),
                left_canon_shape[2].get_center(),
                stroke_width=BOND_HEAVY,
                color=COLOR_ORANGE_RED,
            ).set_z_index(-5)
            
            svd_bond.add_updater(
                lambda mob: mob.put_start_and_end_on(
                    ti_shape[2].get_bottom(),
                    left_canon_shape[2].get_center(),
                )
            )

            left_canon_sweep = AnimationGroup(
                GrowFromEdge(left_canon_shape[0], RIGHT),
                LaggedStart(*[GrowFromEdge(h, RIGHT) for h in left_canon_shape[1]], lag_ratio=0.03),
                GrowFromEdge(left_canon_shape[2], RIGHT),
                Write(left_canon_label_mob),
                lag_ratio=0.0,
            )

            text_swap_anims = []
            if ti_label_mob is not None:
                text_swap_anims.append(Unwrite(ti_label_mob, run_time=move_time / 2))
            text_swap_anims.append(Write(remainder_label_mob, run_time=move_time / 2))
            text_swap = Succession(*text_swap_anims)

            self.play(
                ti_shape.animate.move_to(ti_target),
                text_swap,
                left_canon_sweep,
                Create(svd_bond),
                run_time=move_time,
            )

            if ti_label_mob is not None:
                self.remove(ti_label_mob)
                ti_label_mob.clear_updaters()
            remainder_label_mob.move_to(ti_shape[2].get_center())
            remainder_label_mob.clear_updaters()
            svd_bond.clear_updaters()
            svd_bond.put_start_and_end_on(ti_shape[2].get_bottom(), left_canon_shape[2].get_center())
            
            for line in attached_leg_lines:
                line.clear_updaters()
                line.shift(UP * y_shift)

            ti_shape[0].set_fill("#D3D3D3", opacity=1.0)
            for hatch in ti_shape[1]:
                hatch.set_color(remainder_color)

            self.add(left_canon_shape, left_canon_label_mob, svd_bond)

            return {
                "left_canon_tensor": VGroup(left_canon_shape, left_canon_label_mob),
                "remainder_tensor": VGroup(ti_shape, remainder_label_mob),
                "svd_bond": svd_bond,
                "legs": VGroup(*attached_leg_lines),
                "left_canon_anchors": {
                    "left": left_canon_shape[2].get_left(),
                    "right": left_canon_shape[2].get_right(),
                    "top": left_canon_shape[2].get_top(),
                    "bottom": left_canon_shape[2].get_bottom(),
                    "center": left_canon_shape[2].get_center(),
                },
                "remainder_anchors": {
                    "left": ti_shape[2].get_left(),
                    "right": ti_shape[2].get_right(),
                    "top": ti_shape[2].get_top(),
                    "bottom": ti_shape[2].get_bottom(),
                    "center": ti_shape[2].get_center(),
                },
                "current_wire_y": current_wire_y,
                "neighbor_wire_y": neighbor_wire_y,
            }

        def adjust_circuit_spacing(y_cut, delta, duration=0.5):
            if abs(delta) <= 1e-4:
                return []
            anims = []
            seen_mobs = set()
            stack = list(self.mobjects)
            basic_mobs_to_shift_up = []
            basic_mobs_to_shift_down = []
            lines_to_stretch = []

            while stack:
                mob = stack.pop()
                if id(mob) in seen_mobs: continue
                seen_mobs.add(id(mob))

                if isinstance(mob, Line):
                    s = mob.get_start()
                    e = mob.get_end()
                    is_s_up = s[1] > y_cut
                    is_e_up = e[1] > y_cut
                    
                    if is_s_up and is_e_up:
                        basic_mobs_to_shift_up.append(mob)
                    elif not is_s_up and not is_e_up:
                        basic_mobs_to_shift_down.append(mob)
                    else:
                        lines_to_stretch.append(mob)
                elif getattr(mob, "submobjects", []):
                    stack.extend(mob.submobjects)
                else: 
                    if getattr(mob, "get_center", None) is not None:
                        if mob.get_center()[1] > y_cut:
                            basic_mobs_to_shift_up.append(mob)
                        else:
                            basic_mobs_to_shift_down.append(mob)

            if basic_mobs_to_shift_up:
                group_up = Group(*basic_mobs_to_shift_up)
                anims.append(group_up.animate(run_time=duration).shift(UP * delta / 2))
            if basic_mobs_to_shift_down:
                group_down = Group(*basic_mobs_to_shift_down)
                anims.append(group_down.animate(run_time=duration).shift(DOWN * delta / 2))
            
            for line in lines_to_stretch:
                s = line.get_start()
                e = line.get_end()
                new_s = s + UP * (delta / 2 if s[1] > y_cut else -delta / 2)
                new_e = e + UP * (delta / 2 if e[1] > y_cut else -delta / 2)
                anims.append(line.animate(run_time=duration).put_start_and_end_on(new_s, new_e))
                
            return anims

        def find_all_gates_on_wire(wire_y, tolerance=1e-3):
            gates = []
            if 'copy_items' in locals():
                for pt in copy_items:
                    if hasattr(pt, "gates"):
                        for g in pt.gates:
                            if abs(g.get_center()[1] - wire_y) < tolerance:
                                gates.append(g)
            if 'dt_items' in locals():
                for item in dt_items:
                    if isinstance(item, InvertedPhaseTrain):
                        for g in item.gates:
                            if abs(g.get_center()[1] - wire_y) < tolerance:
                                gates.append(g)
                    elif isinstance(item, HatchedGate):
                        if abs(item.get_center()[1] - wire_y) < tolerance:
                            gates.append(item)
            return sorted(gates, key=lambda g: g.get_center()[0])

        def merge_remainder_upward_with_wire_tensor(
            left_canon_split_info,
            wire_tensor_info,
            output_label=None,
            show_label=False,
            move_time=0.45,
            fuse_time=0.45,
            remove_leg_lines=None,
            extra_leg_lines=None,
        ):
            remainder_shape = left_canon_split_info["remainder_tensor"][0]
            remainder_label = left_canon_split_info["remainder_tensor"][1]
            wire_shape = wire_tensor_info["shape"]
            wire_label = wire_tensor_info.get("label", None)

            remove_leg_lines = list(remove_leg_lines or [])
            extra_leg_lines = list(extra_leg_lines or [])

            left_canon_shape = left_canon_split_info["left_canon_tensor"][0]
            red_bond = left_canon_split_info["svd_bond"]
            leg_lines = list(left_canon_split_info["legs"]) + extra_leg_lines
            
            target_y = wire_shape[2].get_center()[1]
            remainder_target = np.array([remainder_shape[2].get_center()[0], target_y, 0.0])

            # Red bond tracks left-canon (below) and remainder (moving part)
            red_bond.add_updater(
                lambda mob: mob.put_start_and_end_on(
                    remainder_shape[2].get_bottom(),
                    left_canon_shape[2].get_center(),
                )
            )

            for line in leg_lines:
                start_p = line.get_start().copy()
                end_p = line.get_end().copy()
                def rem_up_leg_updater(mob, s=start_p, e=end_p, t=remainder_shape):
                    current_y = t[2].get_center()[1]
                    start_y = left_canon_split_info["remainder_anchors"]["center"][1]
                    dy = current_y - start_y
                    mob.put_start_and_end_on(s + UP * dy, e + UP * dy)
                line.add_updater(rem_up_leg_updater)

            self.play(
                remainder_shape.animate.move_to(remainder_target),
                *[FadeOut(line) for line in remove_leg_lines],
                run_time=move_time
            )

            red_bond.clear_updaters()
            for line in leg_lines: 
                line.clear_updaters()
                dy = target_y - left_canon_split_info["remainder_anchors"]["center"][1]
                line.shift(UP * dy)

            # Final fusion
            final_x_left = min(wire_shape[2].get_left()[0], remainder_shape[2].get_left()[0])
            final_x_right = max(wire_shape[2].get_right()[0], remainder_shape[2].get_right()[0])
            final_wire_y = target_y
            
            ti_color = COLOR_T_US
            final_combined_shape = make_wide_hatched_tensor(final_x_left, final_x_right, final_wire_y, ti_color)
            
            combined_label = None
            if show_label and output_label:
                combined_label = MathTex(output_label, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
                combined_label.move_to(final_combined_shape[2].get_center()).set_z_index(8)

            fuse_anims = [
                ReplacementTransform(wire_shape, final_combined_shape.copy()),
                ReplacementTransform(remainder_shape, final_combined_shape.copy()),
                FadeOut(red_bond),
            ]
            if wire_label: fuse_anims.append(Unwrite(wire_label))
            fuse_anims.append(Unwrite(remainder_label))
            
            if combined_label:
                self.play(*fuse_anims, Write(combined_label), run_time=fuse_time)
            else:
                self.play(*fuse_anims, run_time=fuse_time)

            self.remove(wire_shape, remainder_shape, red_bond)
            if combined_label:
                self.add(final_combined_shape, combined_label)
            else:
                self.add(final_combined_shape)

            return {
                "tensor": VGroup(final_combined_shape, combined_label) if combined_label else VGroup(final_combined_shape),
                "shape": final_combined_shape,
                "label": combined_label,
            }

        def merge_bridge_tensor_with_wire_tensor(
            right_canon_split_info,
            wire_tensor_info,
            output_label=None,
            show_label=False,
            move_time=0.45,
            fuse_time=0.45,
            remove_leg_lines=None,
            extra_leg_lines=None,
        ):
            remainder_shape = right_canon_split_info["remainder_tensor"][0]
            remainder_label = right_canon_split_info["remainder_tensor"][1]
            wire_shape = wire_tensor_info["shape"]
            wire_label = wire_tensor_info.get("label", None)

            remove_leg_lines = list(remove_leg_lines or [])
            extra_leg_lines = list(extra_leg_lines or [])

            right_canon_shape = right_canon_split_info["right_canon_tensor"][0]
            red_bond = right_canon_split_info["svd_bond"]
            raw_leg_lines = list(right_canon_split_info["legs"]) + extra_leg_lines
            dedup_leg_lines = []
            seen_leg_ids = set()
            for line in raw_leg_lines:
                if id(line) not in seen_leg_ids:
                    dedup_leg_lines.append(line)
                    seen_leg_ids.add(id(line))
            leg_lines = [
                line for line in dedup_leg_lines
                if isinstance(line, Line) and line.get_num_points() > 0
            ]
            remove_leg_ids = {
                id(line)
                for line in remove_leg_lines
                if isinstance(line, Line)
            }

            target_y = wire_shape[2].get_center()[1]
            remainder_target = np.array([remainder_shape[2].get_center()[0], target_y, 0.0])

            # Stage A: move remainder down to wire tensor y, with red/white bonds smoothly following.
            red_bond.add_updater(
                lambda mob: mob.put_start_and_end_on(
                    right_canon_shape[2].get_center(),
                    remainder_shape[2].get_top(),
                )
            )

            for line in leg_lines:
                start_point = line.get_start().copy()
                end_x = line.get_end()[0]

                def leg_updater(mob, src=start_point, x_fixed=end_x):
                    if mob.get_num_points() == 0:
                        return
                    src_arr = np.asarray(src)
                    if src_arr.ndim != 1 or src_arr.size != 3:
                        return
                    target = np.array([x_fixed, remainder_shape[2].get_top()[1], 0.0])
                    mob.become(
                        Line(
                            src_arr,
                            target,
                            color=mob.get_color(),
                            stroke_width=mob.stroke_width,
                        ).set_z_index(mob.get_z_index())
                    )

                line.add_updater(leg_updater)

            remainder_label.add_updater(lambda mob: mob.move_to(remainder_shape[2].get_center()))
            self.play(
                Unwrite(remainder_label),
                remainder_shape.animate.move_to(remainder_target),
                run_time=move_time,
            )
            remainder_label.clear_updaters()

            # Stage B: fuse remainder tensor with wire tensor into one bridge tensor.
            x_left = min(remainder_shape[2].get_left()[0], wire_shape[2].get_left()[0]) - 0.02
            x_right = max(remainder_shape[2].get_right()[0], wire_shape[2].get_right()[0]) + 0.02
            y_merge = wire_shape[2].get_center()[1]
            merged_shape = make_wide_hatched_tensor(x_left, x_right, y_merge, COLOR_T_US)

            merged_label = None
            if show_label and output_label is not None:
                merged_label = MathTex(output_label, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
                merged_label.move_to(merged_shape[2].get_center())
                merged_label.set_z_index(8)

            red_bond.clear_updaters()
            red_bond.add_updater(
                lambda mob: mob.put_start_and_end_on(
                    right_canon_shape[2].get_center(),
                    merged_shape[2].get_top(),
                )
            )

            for line in leg_lines:
                line.clear_updaters()
                start_point = line.get_start().copy()
                end_x = line.get_end()[0]

                def leg_merge_updater(mob, src=start_point, x_fixed=end_x):
                    if mob.get_num_points() == 0:
                        return
                    src_arr = np.asarray(src)
                    if src_arr.ndim != 1 or src_arr.size != 3:
                        return
                    target = np.array([x_fixed, merged_shape[2].get_top()[1], 0.0])
                    mob.become(
                        Line(
                            src_arr,
                            target,
                            color=mob.get_color(),
                            stroke_width=mob.stroke_width,
                        ).set_z_index(mob.get_z_index())
                    )

                line.add_updater(leg_merge_updater)

            remainder_target_copy = merged_shape.copy()
            wire_target_copy = merged_shape.copy()

            transform_anims = [
                ReplacementTransform(remainder_shape, remainder_target_copy),
                ReplacementTransform(wire_shape, wire_target_copy),
            ]
            if wire_label is not None:
                transform_anims.append(Unwrite(wire_label))
            if merged_label is not None:
                transform_anims.append(Write(merged_label))
            for line in leg_lines:
                if id(line) in remove_leg_ids:
                    transform_anims.append(FadeOut(line))

            self.play(*transform_anims, run_time=fuse_time)

            self.remove(remainder_target_copy, wire_target_copy)

            self.add(merged_shape)
            if merged_label is not None:
                self.add(merged_label)

            red_bond.clear_updaters()
            red_bond.put_start_and_end_on(right_canon_shape[2].get_center(), merged_shape[2].get_top())
            for line in leg_lines:
                line.clear_updaters()
                if id(line) in remove_leg_ids:
                    continue
                start_point = line.get_start().copy()
                end_x = line.get_end()[0]
                start_arr = np.asarray(start_point)
                if start_arr.ndim == 1 and start_arr.size == 3 and line.get_num_points() > 0:
                    line.become(
                        Line(
                            start_arr,
                            np.array([end_x, merged_shape[2].get_top()[1], 0.0]),
                            color=line.get_color(),
                            stroke_width=line.stroke_width,
                        ).set_z_index(line.get_z_index())
                    )

            # Add red_bond back to scene so it can be found by compress_oc_recursive
            self.add(red_bond)

            remaining_leg_lines = [
                line for line in leg_lines
                if id(line) not in remove_leg_ids
            ]

            return {
                "tensor": VGroup(merged_shape, merged_label) if merged_label is not None else VGroup(merged_shape),
                "shape": merged_shape,
                "label": merged_label,
                "svd_bond": red_bond,
                "legs": VGroup(*remaining_leg_lines),
                "anchors": {
                    "left": merged_shape[2].get_left(),
                    "right": merged_shape[2].get_right(),
                    "top": merged_shape[2].get_top(),
                    "bottom": merged_shape[2].get_bottom(),
                    "center": merged_shape[2].get_center(),
                },
            }

        def finalize_single_leg_tensor_step(
            hd_gate,
            control_dot,
            control_line,
            target_tensor_info,
            pivot_x,
            unitary_label="U_{3}",
            color=COLOR_U,
            move_time=0.5,
            reshape_time=0.5,
        ):
            target_shape = target_tensor_info["shape"]
            target_label = target_tensor_info.get("label", None)

            # Stage A: H_d moves to control, morphs hatch/style, and swaps label to U_3 during motion.
            hd_body = VGroup(hd_gate.bg, hd_gate.hatching, hd_gate.box)
            control_center = control_dot.get_center()
            u3_gate = HatchedGate(unitary_label, control_center[0], control_center[1], color)
            u3_body = VGroup(u3_gate.bg, u3_gate.hatching, u3_gate.box)
            u3_label = MathTex(unitary_label, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
            u3_label.set_z_index(8)

            hd_gate.label.clear_updaters()
            hd_gate.label.move_to(hd_body[2].get_center())
            u3_label.move_to(control_center)

            self.play(
                hd_body.animate.move_to(control_center),
                Transform(hd_body, u3_body),
                FadeOut(control_dot),
                ReplacementTransform(hd_gate.label, u3_label),
                run_time=0.35,
            )

            self.remove(hd_gate)
            self.add(hd_body, u3_label)
            u3_box = hd_body[2]
            u3_label.add_updater(lambda mob: mob.move_to(u3_box.get_center()))

            # Keep control bond attached to U_3 and T_2 during motion.
            control_line.add_updater(
                lambda mob: mob.put_start_and_end_on(
                    u3_box.get_bottom(),
                    np.array([
                        np.clip(
                            u3_box.get_center()[0],
                            target_shape[2].get_left()[0] + 0.02,
                            target_shape[2].get_right()[0] - 0.02,
                        ),
                        target_shape[2].get_top()[1],
                        0.0,
                    ]),
                )
            )

            # Stage B/C: move U_3 to V_1 pivot axis while reshaping wide T_2 -> square T_2.
            u3_target = np.array([pivot_x, control_center[1], 0.0])
            t2_center = target_shape[2].get_center()
            square_t2 = make_wide_hatched_tensor(
                t2_center[0] - 0.5 * GATE_WIDTH,
                t2_center[0] + 0.5 * GATE_WIDTH,
                t2_center[1],
                COLOR_T_US,
            )

            if target_label is not None:
                target_label.add_updater(lambda mob: mob.move_to(target_shape[2].get_center()))

            self.play(
                AnimationGroup(
                    hd_body.animate(run_time=move_time).move_to(u3_target),
                    Transform(target_shape, square_t2, run_time=reshape_time),
                    lag_ratio=0.0,
                )
            )

            control_line.clear_updaters()
            control_line.put_start_and_end_on(
                u3_box.get_bottom(),
                np.array([
                    np.clip(
                        u3_box.get_center()[0],
                        target_shape[2].get_left()[0] + 0.02,
                        target_shape[2].get_right()[0] - 0.02,
                    ),
                    target_shape[2].get_top()[1],
                    0.0,
                ]),
            )
            u3_label.clear_updaters()
            u3_label.move_to(u3_box.get_center())

            if target_label is not None:
                target_label.clear_updaters()
                target_label.move_to(target_shape[2].get_center())

            return {
                "u3_tensor": VGroup(hd_body, u3_label),
                "tensor": VGroup(target_shape, target_label) if target_label is not None else VGroup(target_shape),
                "shape": target_shape,
                "label": target_label,
                "control_line": control_line,
                "anchors": {
                    "left": target_shape[2].get_left(),
                "right": target_shape[2].get_right(),
                    "top": target_shape[2].get_top(),
                    "bottom": target_shape[2].get_bottom(),
                    "center": target_shape[2].get_center(),
                },
            }

        def compress_oc_recursive(
            train_idx,
            t_current,
            canon_map,
            direction="down",
            previous_canon=None,
            blue_bond_to_current_t=None,
            show_pulse=True,
            y_override_list=None,
            label_suffix="",
            prime_canon_indices=None,
        ):
            def _shape_box(shape_mob):
                if isinstance(shape_mob, VGroup) and len(shape_mob) >= 3:
                    return shape_mob[2]
                if isinstance(shape_mob, VGroup) and len(shape_mob) >= 1:
                    return shape_mob[0]
                return shape_mob

            def _iter_scene_lines(root_mobjects):
                stack = list(root_mobjects)
                seen = set()
                while stack:
                    mob = stack.pop()
                    mid = id(mob)
                    if mid in seen:
                        continue
                    seen.add(mid)
                    if isinstance(mob, Line):
                        yield mob
                    submobs = getattr(mob, "submobjects", [])
                    if submobs:
                        stack.extend(submobs)

            def _find_inter_tensor_bond_candidates(box_a, box_b, tol=0.5):
                anchor_a = _boundary_anchor_for_point(box_a, box_b.get_center())
                anchor_b = _boundary_anchor_for_point(box_b, box_a.get_center())
                candidates = []

                for line in _iter_scene_lines(self.mobjects):
                    if line.get_num_points() == 0:
                        continue
                    start = line.get_start()
                    end = line.get_end()
                    score_direct = np.linalg.norm(start - anchor_a) + np.linalg.norm(end - anchor_b)
                    score_swap = np.linalg.norm(start - anchor_b) + np.linalg.norm(end - anchor_a)
                    score = min(score_direct, score_swap)
                    if score <= tol:
                        candidates.append((score, line))

                candidates.sort(key=lambda item: item[0])
                return [line for _, line in candidates]

            if direction not in ("down", "up"):
                raise ValueError("direction must be either 'down' or 'up'")

            t_current_box = _shape_box(t_current["shape"])
            axis_x = t_current_box.get_center()[0]
            oc_bond_width = BOND_LIGHT
            next_idx = train_idx + 1 if direction == "down" else train_idx - 1

            if next_idx < 1 or next_idx > n_main:
                return t_current, previous_canon, blue_bond_to_current_t

            canon_next = canon_map.get(next_idx)
            if canon_next is None:
                return t_current, previous_canon, blue_bond_to_current_t

            if show_pulse:
                pulse_highlight_box(VGroup(t_current["shape"], canon_next["shape"]))

            canon_next_box = _shape_box(canon_next["shape"])
            y_top = max(t_current_box.get_top()[1], canon_next_box.get_top()[1])
            y_bottom = min(t_current_box.get_bottom()[1], canon_next_box.get_bottom()[1])

            effective_y = y_override_list if y_override_list is not None else [main_y(i) for i in range(n_main)]

            merged_pair = make_tall_hatched_tensor(axis_x, y_top, y_bottom, COLOR_T_US)
            merge_target_t = merged_pair.copy()
            merge_target_canon = merged_pair.copy()

            merge_anims = [
                ReplacementTransform(t_current["shape"], merge_target_t),
                ReplacementTransform(canon_next["shape"], merge_target_canon),
            ]

            explicit_local_bond = canon_next.get("svd_bond")
            if explicit_local_bond is None:
                explicit_local_bond = t_current.get("svd_bond")
            if (
                isinstance(explicit_local_bond, Line)
                and explicit_local_bond.get_num_points() > 0
            ):
                explicit_local_bond.clear_updaters()
                explicit_local_bond.put_start_and_end_on(
                    _boundary_anchor_for_point(canon_next_box, t_current_box.get_center()),
                    _boundary_anchor_for_point(t_current_box, canon_next_box.get_center()),
                )
                self.add(explicit_local_bond)
            else:
                explicit_local_bond = None

            local_bond_candidates = _find_inter_tensor_bond_candidates(t_current_box, canon_next_box)
            local_bond = explicit_local_bond if explicit_local_bond is not None else (local_bond_candidates[0] if local_bond_candidates else None)
            if local_bond is not None:
                local_bond.clear_updaters()
            for extra_bond in local_bond_candidates:
                if local_bond is not None and id(extra_bond) == id(local_bond):
                    continue
                extra_bond.clear_updaters()
                merge_anims.append(FadeOut(extra_bond))
            if local_bond is not None:
                local_bond_target = Line(
                    local_bond.get_start(),
                    local_bond.get_end(),
                    color=COLOR_LIGHT_BLUE,
                    stroke_width=oc_bond_width,
                ).set_z_index(local_bond.get_z_index())
                merge_anims.append(Transform(local_bond, local_bond_target))

            if t_current.get("label") is not None:
                merge_anims.append(Unwrite(t_current["label"]))
            if canon_next.get("label") is not None:
                merge_anims.append(Unwrite(canon_next["label"]))

            blue_bond_to_merged = None
            if previous_canon is not None and blue_bond_to_current_t is not None:
                prev_box = _shape_box(previous_canon["shape"])
                merged_box = merged_pair[2]
                merged_anchor = _boundary_anchor_for_point(merged_box, prev_box.get_center())
                prev_anchor = _boundary_anchor_for_point(prev_box, merged_box.get_center())
                blue_bond_to_merged = Line(prev_anchor, merged_anchor, color=COLOR_LIGHT_BLUE, stroke_width=oc_bond_width).set_z_index(0)
                merge_anims.append(ReplacementTransform(blue_bond_to_current_t, blue_bond_to_merged))

            self.play(*merge_anims, run_time=0.42)
            self.remove(merge_target_t, merge_target_canon)
            self.add(merged_pair)
            self.wait(0.15)

            if direction == "down":
                canon_idx = train_idx
                t_next_idx = train_idx + 1
                current_suffix = (
                    label_suffix
                    if prime_canon_indices is not None and canon_idx in prime_canon_indices
                    else ""
                )
                canon_shape = make_square_hatched_tensor(axis_x, effective_y[canon_idx - 1], COLOR_U)
                t_next_shape = make_square_hatched_tensor(axis_x, effective_y[t_next_idx - 1], COLOR_T_US)
                canon_label = MathTex(f"U_{{{canon_idx}}}{{{current_suffix}}}", color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
                t_next_label = MathTex(f"T_{{{t_next_idx}}}", color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
            else:
                canon_idx = train_idx
                t_next_idx = train_idx - 1
                current_suffix = (
                    label_suffix
                    if prime_canon_indices is not None and canon_idx in prime_canon_indices
                    else ""
                )
                canon_shape = make_square_hatched_tensor(axis_x, effective_y[canon_idx - 1], COLOR_PURPLE)
                t_next_shape = make_square_hatched_tensor(axis_x, effective_y[t_next_idx - 1], COLOR_T_US)
                canon_label = MathTex(f"V_{{{canon_idx}}}{{{current_suffix}}}", color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
                t_next_label = MathTex(f"T_{{{t_next_idx}}}", color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)

            canon_label.move_to(canon_shape[2].get_center())
            t_next_label.move_to(t_next_shape[2].get_center())
            canon_label.set_z_index(8)
            t_next_label.set_z_index(8)

            canon_box = canon_shape[2]
            t_box = t_next_shape[2]

            if canon_box.get_center()[1] >= t_box.get_center()[1]:
                bond_start, bond_end = canon_box.get_bottom(), t_box.get_top()
            else:
                bond_start, bond_end = canon_box.get_top(), t_box.get_bottom()

            vertical_blue_bond = Line(bond_start, bond_end, color=COLOR_LIGHT_BLUE, stroke_width=oc_bond_width).set_z_index(0)

            split_anims = [
                ReplacementTransform(merged_pair, canon_shape),
                Create(t_next_shape),
                Write(canon_label),
                Write(t_next_label),
            ]
            if local_bond is not None:
                local_bond.clear_updaters()
                split_anims.append(ReplacementTransform(local_bond, vertical_blue_bond))
            else:
                split_anims.append(Create(vertical_blue_bond))

            new_upstream_bond = None
            if previous_canon is not None and blue_bond_to_merged is not None:
                prev_box = _shape_box(previous_canon["shape"])
                prev_anchor = _boundary_anchor_for_point(prev_box, canon_box.get_center())
                canon_anchor = _boundary_anchor_for_point(canon_box, prev_box.get_center())
                new_upstream_bond = Line(prev_anchor, canon_anchor, color=COLOR_LIGHT_BLUE, stroke_width=oc_bond_width).set_z_index(0)
                split_anims.append(ReplacementTransform(blue_bond_to_merged, new_upstream_bond))

            self.play(*split_anims, run_time=0.75)
            self.add(canon_shape, canon_label, t_next_shape, t_next_label, vertical_blue_bond)
            if new_upstream_bond is not None:
                self.add(new_upstream_bond)

            canon_map.pop(next_idx, None)

            t_next_info = {
                "shape": t_next_shape,
                "label": t_next_label,
                "wire_idx": t_next_idx,
            }
            canon_info = {
                "shape": canon_shape,
                "label": canon_label,
                "wire_idx": canon_idx,
                "upstream_canon": previous_canon,
                "upstream_bond": new_upstream_bond,
            }

            return compress_oc_recursive(
                train_idx=t_next_idx,
                t_current=t_next_info,
                canon_map=canon_map,
                direction=direction,
                previous_canon=canon_info,
                blue_bond_to_current_t=vertical_blue_bond,
                show_pulse=show_pulse,
                y_override_list=y_override_list,
                label_suffix=label_suffix,
                prime_canon_indices=prime_canon_indices,
            )

        def set_objects_opacity(objects, opacity, run_time=0.6):
            flat_objects = []
            for obj in objects:
                if isinstance(obj, (list, tuple, VGroup, Group)):
                    for inner in obj:
                        flat_objects.append(inner)
                else:
                    flat_objects.append(obj)
            if flat_objects:
                if opacity >= 0.999:
                    restore_anims = []
                    fallback_anims = []
                    for obj in flat_objects:
                        if getattr(obj, "_opacity_state_saved", False):
                            restore_anims.append(Restore(obj))
                            obj._opacity_state_saved = False
                        else:
                            fallback_anims.append(obj.animate.set_opacity(opacity))
                    self.play(*(restore_anims + fallback_anims), run_time=run_time)
                else:
                    dim_anims = []
                    for obj in flat_objects:
                        if not getattr(obj, "_opacity_state_saved", False):
                            obj.save_state()
                            obj._opacity_state_saved = True
                        dim_anims.append(obj.animate.set_opacity(opacity))
                    self.play(*dim_anims, run_time=run_time)

        def pulse_highlight_box(target_mob, pulse_count=2, pulse_time=0.25, low_alpha=0.25, high_alpha=0.9):
            highlight = RoundedRectangle(
                corner_radius=0.12,
                width=target_mob.width + 0.25,
                height=target_mob.height + 0.25,
                color=RED, stroke_width=3,
                fill_color=RED, fill_opacity=0.075,
            ).move_to(target_mob.get_center()).set_z_index(20)
            highlight.set_stroke(opacity=low_alpha)
            self.add(highlight)
            for _ in range(pulse_count):
                self.play(highlight.animate.set_stroke(opacity=high_alpha), run_time=pulse_time)
                self.play(highlight.animate.set_stroke(opacity=low_alpha), run_time=pulse_time)
            self.remove(highlight)

        # ================================================================
        # PHASE 0: CIRCUIT LAYOUT AND REVEAL
        # ================================================================
        DT_SLOTS, COPY_SLOTS = 7, 3
        SECTION_GAP = 0.5
        WIRE_EDGE_MARGIN = MIN_GATE_GAP
        total_width = (DT_SLOTS + COPY_SLOTS) * GATE_CENTER_SPACING + SECTION_GAP + 2 * WIRE_EDGE_MARGIN
        wire_start_x = CIRCUIT_CENTER_X - total_width / 2

        first_x = wire_start_x + WIRE_EDGE_MARGIN + GATE_WIDTH / 2
        dt_x = [first_x + i * GATE_CENTER_SPACING for i in range(DT_SLOTS)]
        copy_x = [dt_x[-1] + GATE_CENTER_SPACING + SECTION_GAP + i * GATE_CENTER_SPACING for i in range(COPY_SLOTS)]

        def main_y(i): return wire_y_positions[2 * i]
        def copy_y(i): return wire_y_positions[2 * i + 1]

        wires = []
        for wi in range(n_total):
            w = QuantumWire(wire_y_positions[wi], total_width, wire_start_x, wi, is_copy=(wi % 2 == 1))
            wires.append(w)
            self.play(Create(w.line), FadeIn(w.start_circ), FadeIn(w.end_circ), Write(w.label), run_time=0.15)

        dt_items = []
        for i in range(n_main):
            hd = HatchedGate("H_d", dt_x[2*i], main_y(i), COLOR_HD)
            dt_items.append(hd)
            if i < n_main - 1:
                targets = [(main_y(tj), f"R_{{{tj+1}{i+2}}}", [COLOR_BROWN_1, COLOR_BROWN_2, COLOR_BROWN_3][min(i-tj, 2)], True) for tj in range(i + 1)]
                dt_items.append(InvertedPhaseTrain(dt_x[2*i+1], main_y(i+1), targets))

        copy_items = []
        for j in range(n_main - 1):
            targets = [(main_y(ti), f"R_{{{ti+1}{j+1}}}", [COLOR_BROWN_1, COLOR_BROWN_2, COLOR_BROWN_3][min(ti-j-1, 2)], True) for ti in range(j + 1, n_main)]
            copy_items.append(PhaseTrain(copy_x[j], copy_y(j), targets))



        def reveal(item):
            if isinstance(item, (QuantumGate, HatchedGate)):
                anims = []
                if hasattr(item, 'bg'):
                    item.bg.save_state()
                    item.bg.stretch_to_fit_height(0.01)
                    anims.append(Restore(item.bg))

                if hasattr(item, 'hatching'):
                    for hl in item.hatching:
                        hl.save_state()
                        hl.stretch_to_fit_height(0.01)
                        anims.append(Restore(hl))
                
                item.box.save_state()
                item.box.stretch_to_fit_height(0.01)
                anims.append(Restore(item.box))
                anims.append(Write(item.label))
                self.play(*anims, run_time=0.1)
            else:
                item.remove(item.cline)
                self.play(GrowFromCenter(item.dot), Create(item.cline), run_time=0.1)
                for g in item.gates:
                    reveal(g)
                item.add(item.cline)

        # Draw the full DT circuit (DT ∘ Copy) before any compression phase.
        for item in [*dt_items, *copy_items]:
            reveal(item)

        self.wait(1.0)

        # Focus stage: keep only main register + DT section prominent.
        copy_wire_parts = []
        for w in wires:
            if w.is_copy:
                copy_wire_parts.extend([w.line, w.start_circ, w.end_circ, w.label])
        set_objects_opacity([copy_wire_parts, copy_items], opacity=0.18, run_time=0.8)

        # ================================================================
        # PHASE 1 ~ PHASE 3: DT MPO BUILD AND PRE-COPY CANONICALIZATION
        # ================================================================
        # Pilot TT merge on n1 using the first two DT TT blocks: H_d + R_{12} + R_{13}
        dt_bond_anchor_registry = {}
        n1_dt_targets = get_dt_wire_target_gates(dt_items, main_y(0), train_indices=[1, 3])
        if len(n1_dt_targets) != 2:
            raise ValueError("Expected exactly R_{12} and R_{13} on n_1 for the pilot TT merge")

        pulse_highlight_box(VGroup(dt_items[0], *n1_dt_targets), pulse_count=2, pulse_time=0.18)

        n1_tt_merge_result = merge_wire_tt_block(
            wire_y=main_y(0),
            hd_gate=dt_items[0],
            target_gates=n1_dt_targets,
            ti_label="T_{1}",
            ti_color=COLOR_T_US,
            bond_anchor_registry=dt_bond_anchor_registry,
            show_label=False,
        )

        n1_split_result = split_tensor_to_right_canon(
            tensor_info=n1_tt_merge_result,
            current_wire_y=main_y(0),
            neighbor_wire_y=main_y(1),
            right_canon_label="Q_{1}",
            remainder_label="R_{1}",
            attached_leg_lines=[
                dt_items[1].cline,
                dt_items[3].cline,
            ],
            right_canon_color=COLOR_U,
            remainder_color=COLOR_T_US,
            move_time=0.75,
        )

        # Next step: merge n2 control + H_d + R_{23} into n2 block
        n2_r23_targets = get_dt_wire_target_gates(dt_items, main_y(1), train_indices=[3])
        if len(n2_r23_targets) != 1:
            raise ValueError("Expected exactly one R_{23} on n_2 for the n_2 merge step")

        n2_merge_sources = [dt_items[1].dot, n2_r23_targets[0]]
        pulse_highlight_box(VGroup(dt_items[2], *n2_merge_sources, n1_split_result["remainder_tensor"][0]), pulse_count=2, pulse_time=0.16)

        n2_tt_merge_result = merge_wire_tt_block(
            wire_y=main_y(1),
            hd_gate=dt_items[2],
            target_gates=n2_merge_sources,
            ti_label="T_{2}",
            ti_color=COLOR_T_US,
            bond_anchor_registry=dt_bond_anchor_registry,
            show_label=False,
        )

        # Bridge-merge step: merge R_1 with the newly combined n2 tensor → R_2
        n2_bridge_merge_result = merge_bridge_tensor_with_wire_tensor(
            right_canon_split_info=n1_split_result,
            wire_tensor_info=n2_tt_merge_result,
            output_label="R_{2}",
            show_label=True,
            move_time=0.42,
            fuse_time=0.45,
            remove_leg_lines=[dt_items[1].cline],
        )

        # Finalize at n3 using control (do NOT merge H_d with R_{34}) → R_3
        pulse_highlight_box(VGroup(dt_items[4], dt_items[3].dot), pulse_count=2, pulse_time=0.16)
        n3_finalize_result = finalize_single_leg_tensor_step(
            hd_gate=dt_items[4],
            control_dot=dt_items[3].dot,
            control_line=dt_items[3].cline,
            target_tensor_info=n2_bridge_merge_result,
            pivot_x=n1_split_result["right_canon_anchors"]["center"][0],
            unitary_label="R_{3}",
            color=COLOR_T_US,
            move_time=0.48,
            reshape_time=0.45,
        )

        # Split R_2 → Q_2 (stays at n2, left-canon color) + R_2' (moves toward n3)
        r2_tensor_info = {
            "shape": n3_finalize_result["shape"],
            "label": n3_finalize_result["label"],
        }
        r2_split_result = split_tensor_to_right_canon(
            tensor_info=r2_tensor_info,
            current_wire_y=main_y(1),
            neighbor_wire_y=main_y(2),
            right_canon_label="Q_{2}",
            remainder_label="R_{2}'",
            right_canon_color=COLOR_U,
            remainder_color=COLOR_T_US,
            attached_leg_lines=[n3_finalize_result["control_line"]],
            move_time=0.72,
        )

        # Merge R_2' + R_3 → T_3 (orthogonality center)
        r3_body = n3_finalize_result["u3_tensor"][0]
        r3_label_mob = n3_finalize_result["u3_tensor"][1]
        r3_tensor_info = {"shape": r3_body, "label": r3_label_mob}
        t3_merge_result = merge_bridge_tensor_with_wire_tensor(
            right_canon_split_info=r2_split_result,
            wire_tensor_info=r3_tensor_info,
            output_label="T_{3}",
            show_label=True,
            move_time=0.42,
            fuse_time=0.45,
            remove_leg_lines=[n3_finalize_result["control_line"]],
        )

        # OC recursive sweep upward from T_3 with Q_1, Q_2 → produces V_3, V_2, T_1
        oc_t_start_pre = {
            "shape": t3_merge_result["shape"],
            "label": t3_merge_result["label"],
            "wire_idx": 3,
        }
        oc_canon_map_pre = {
            1: {
                "shape": n1_split_result["right_canon_tensor"][0],
                "label": n1_split_result["right_canon_tensor"][1],
                "wire_idx": 1,
            },
            2: {
                "shape": r2_split_result["right_canon_tensor"][0],
                "label": r2_split_result["right_canon_tensor"][1],
                "wire_idx": 2,
            },
        }
        t1_after_pre, v2_after_pre, t1_v2_bond = compress_oc_recursive(
            train_idx=3,
            t_current=oc_t_start_pre,
            canon_map=oc_canon_map_pre,
            direction="up",
        )

        self.wait(0.6)

        # PHASE 3: Continue sequence at n=4
        v3_after_pre = v2_after_pre["upstream_canon"]
        v2_v3_bond = v2_after_pre["upstream_bond"]
        
        dt_train4 = dt_items[5]
        r14 = dt_train4.gates[0]
        r24 = dt_train4.gates[1]
        r34 = dt_train4.gates[2]
        control4 = dt_train4.dot
        control_line4 = dt_train4.cline
        hd4 = dt_items[6]

        # 1. Move T1V2V3 closer to R14 R24 R34 control
        group_to_move = VGroup(
            t1_after_pre["shape"], t1_after_pre["label"],
            v2_after_pre["shape"], v2_after_pre["label"],
            v3_after_pre["shape"], v3_after_pre["label"],
            t1_v2_bond, v2_v3_bond
        )
        target_x = r14.get_center()[0] - 0.75
        shift_amount = target_x - group_to_move.get_center()[0]
        self.play(group_to_move.animate.shift(RIGHT * shift_amount), run_time=0.8)

        # Keep labels attached to shapes so merge_wire_tt_block can handle them automatically
        t1_after_pre["shape"].label = t1_after_pre["label"]
        v2_after_pre["shape"].label = v2_after_pre["label"]
        v3_after_pre["shape"].label = v3_after_pre["label"]

        # 2. T1 merges with R14 -> Q1, R1
        pulse_highlight_box(VGroup(t1_after_pre["shape"], r14), pulse_count=2, pulse_time=0.16)
        n4_w0_merge = merge_wire_tt_block(
            wire_y=main_y(0),
            hd_gate=None,
            target_gates=[t1_after_pre["shape"], r14],
            ti_label="temp",
            ti_color=COLOR_T_US,
            show_label=False,
            attached_leg_lines=[t1_v2_bond, control_line4]
        )
        n4_r1_split = split_tensor_to_right_canon(
            tensor_info=n4_w0_merge,
            current_wire_y=main_y(0),
            neighbor_wire_y=main_y(1),
            right_canon_label="Q_{1}",
            remainder_label="R_{1}",
            right_canon_color=COLOR_U,
            remainder_color=COLOR_T_US,
            attached_leg_lines=[t1_v2_bond, control_line4],
            move_time=0.75
        )

        # 3. V2 combines with R24 -> R1 merges to give R2
        pulse_highlight_box(VGroup(v2_after_pre["shape"], r24, n4_r1_split["remainder_tensor"][0]), pulse_count=2, pulse_time=0.16)
        n4_w1_merge = merge_wire_tt_block(
            wire_y=main_y(1),
            hd_gate=None,
            target_gates=[v2_after_pre["shape"], r24],
            ti_label="temp",
            ti_color=COLOR_T_US,
            show_label=False,
            attached_leg_lines=[v2_v3_bond]
        )
        n4_r2_bridge = merge_bridge_tensor_with_wire_tensor(
            right_canon_split_info=n4_r1_split,
            wire_tensor_info=n4_w1_merge,
            output_label="R_{2}",
            show_label=True,
            move_time=0.42,
            fuse_time=0.45,
            remove_leg_lines=[t1_v2_bond],
            extra_leg_lines=[control_line4]
        )

        # 4. R2 splits -> Q2, R2'. V3 merges with R34. R2' bridges V3+R34 -> R3
        n4_r2_split = split_tensor_to_right_canon(
            tensor_info=n4_r2_bridge,
            current_wire_y=main_y(1),
            neighbor_wire_y=main_y(2),
            right_canon_label="Q_{2}",
            remainder_label="R_{2}'",
            right_canon_color=COLOR_U,
            remainder_color=COLOR_T_US,
            attached_leg_lines=[v2_v3_bond, control_line4],
            move_time=0.72
        )
        pulse_highlight_box(VGroup(v3_after_pre["shape"], r34, n4_r2_split["remainder_tensor"][0]), pulse_count=2, pulse_time=0.16)
        n4_w2_merge = merge_wire_tt_block(
            wire_y=main_y(2),
            hd_gate=None,
            target_gates=[v3_after_pre["shape"], r34],
            ti_label="temp",
            ti_color=COLOR_T_US,
            show_label=False,
            attached_leg_lines=[]
        )
        n4_r3_bridge = merge_bridge_tensor_with_wire_tensor(
            right_canon_split_info=n4_r2_split,
            wire_tensor_info=n4_w2_merge,
            output_label="R_{3}",
            show_label=True,
            move_time=0.42,
            fuse_time=0.45,
            remove_leg_lines=[v2_v3_bond],
            extra_leg_lines=[control_line4]
        )

        # 5. Finalize Hd with control -> R4
        pulse_highlight_box(VGroup(hd4, control4), pulse_count=2, pulse_time=0.16)
        n4_r4_finalize = finalize_single_leg_tensor_step(
            hd_gate=hd4,
            control_dot=control4,
            control_line=control_line4,
            target_tensor_info=n4_r3_bridge,
            pivot_x=n4_r2_split["right_canon_anchors"]["center"][0],
            unitary_label="R_{4}",
            color=COLOR_T_US,
            move_time=0.48,
            reshape_time=0.45
        )

        # 6. QR of R3 to Q3, R4_prime and merge with R4 -> T4
        r3_tensor_info = {
            "shape": n4_r4_finalize["shape"],
            "label": n4_r4_finalize["label"]
        }
        n4_r3_split = split_tensor_to_right_canon(
            tensor_info=r3_tensor_info,
            current_wire_y=main_y(2),
            neighbor_wire_y=main_y(3),
            right_canon_label="Q_{3}",
            remainder_label="R_{4}'",
            right_canon_color=COLOR_U,
            remainder_color=COLOR_T_US,
            attached_leg_lines=[n4_r4_finalize["control_line"]],
            move_time=0.72
        )
        
        r4_finalize_info = {
            "shape": n4_r4_finalize["u3_tensor"][0],
            "label": n4_r4_finalize["u3_tensor"][1]
        }
        n4_t4_merge = merge_bridge_tensor_with_wire_tensor(
            right_canon_split_info=n4_r3_split,
            wire_tensor_info=r4_finalize_info,
            output_label="T_{4}",
            show_label=True,
            move_time=0.42,
            fuse_time=0.45,
            remove_leg_lines=[n4_r4_finalize["control_line"]]
        )

        # 7. Compression sweep to the top
        oc_t_start_final = {
            "shape": n4_t4_merge["shape"],
            "label": n4_t4_merge["label"],
            "wire_idx": 4, 
        }
        oc_canon_map_final = {
            1: {
                "shape": n4_r1_split["right_canon_tensor"][0],
                "label": n4_r1_split["right_canon_tensor"][1],
                "wire_idx": 1,
            },
            2: {
                "shape": n4_r2_split["right_canon_tensor"][0],
                "label": n4_r2_split["right_canon_tensor"][1],
                "wire_idx": 2,
            },
            3: {
                "shape": n4_r3_split["right_canon_tensor"][0],
                "label": n4_r3_split["right_canon_tensor"][1],
                "wire_idx": 3,
            },
        }
        final_t, final_v2, final_bond = compress_oc_recursive(
            train_idx=4,
            t_current=oc_t_start_final,
            canon_map=oc_canon_map_final,
            direction="up"
        )

        self.wait(0.6)
        
        # Fix wire z-index: push all wires to the very back
        for w in wires:
            w.line.set_z_index(-10)
            w.start_circ.set_z_index(-10)
            w.end_circ.set_z_index(-10)
            w.label.set_z_index(-9)

        # ================================================================
        # PHASE 4: DOWNWARD SWEEP INTO THE COPY BLOCK
        # ================================================================

        # Use the actual Phase 3 output as Phase 4 input (no fast-forward reconstruction).
        t1_shape = final_t["shape"]
        t1_label = final_t["label"]
        v2_shape = final_v2["shape"]
        v2_label = final_v2["label"]
        bond12 = final_bond
        v3_info = final_v2["upstream_canon"]
        bond23 = final_v2["upstream_bond"]
        v4_info = v3_info["upstream_canon"] if isinstance(v3_info, dict) else None
        bond34 = v3_info["upstream_bond"] if isinstance(v3_info, dict) else None
        v3_shape = v3_info["shape"] if isinstance(v3_info, dict) else None
        v3_label = v3_info["label"] if isinstance(v3_info, dict) else None
        v4_shape = v4_info["shape"] if isinstance(v4_info, dict) else None
        v4_label = v4_info["label"] if isinstance(v4_info, dict) else None

        group_members = [
            t1_shape, t1_label,
            v2_shape, v2_label,
            v3_shape, v3_label,
            v4_shape, v4_label,
            bond12, bond23, bond34,
        ]
        group_to_move = VGroup(*[mob for mob in group_members if mob is not None])

        # 1. Make the copy block visible again
        set_objects_opacity([copy_wire_parts, copy_items], opacity=1.0, run_time=0.6)

        # 2. Move tensor closer to copy block
        first_copy_pt = copy_items[0]
        target_copy_x = first_copy_pt.dot.get_center()[0] - 0.75
        shift_amount_copy = target_copy_x - group_to_move.get_center()[0]
        self.play(group_to_move.animate.shift(RIGHT * shift_amount_copy), run_time=0.8)

        # 3. Downward sweep: T1 down to T4 with U1 U2 U3 on top
        oc_t_start_down = {
            "shape": t1_shape,
            "label": t1_label,
            "wire_idx": 1,
        }
        oc_canon_map_down = {
            2: {"shape": v2_shape, "label": v2_label, "wire_idx": 2, "svd_bond": bond12},
            3: {"shape": v3_shape, "label": v3_label, "wire_idx": 3, "svd_bond": bond23},
            4: {"shape": v4_shape, "label": v4_label, "wire_idx": 4, "svd_bond": bond34},
        }

        # The function `compress_oc_recursive` natively finds bonds!
        final_t_down, u_info_down, last_bond_down = compress_oc_recursive(
            train_idx=1,
            t_current=oc_t_start_down,
            canon_map=oc_canon_map_down,
            direction="down",
            show_pulse=False
        )

        self.wait(0.6)

        # ================================================================
        # PHASE 5: UPWARD QR OVER FIRST COPY BLOCK
        # ================================================================

        # 1. Access components from the DOWNWARD sweep
        # final_t_down is T4
        # u3_down is the previous_canon from the final call
        u3_down = u_info_down
        bond34_down = last_bond_down
        
        u2_down = u3_down["upstream_canon"]
        bond23_down = u3_down["upstream_bond"]
        
        u1_down = u2_down["upstream_canon"]
        bond12_down = u2_down["upstream_bond"]
        
        t4_down = final_t_down
        # Attach labels to shapes so merge_wire_tt_block handles their Unwrite
        t4_down["shape"].label = t4_down["label"]
        u3_down["shape"].label = u3_down["label"]
        u2_down["shape"].label = u2_down["label"]
        
        # Gates on the FIRST COPY TRAIN (copy_items[0])
        first_copy_train = copy_items[0]
        r21 = first_copy_train.gates[0]
        r31 = first_copy_train.gates[1]
        r41 = first_copy_train.gates[2]
        copy_dot1 = first_copy_train.dot
        copy_cline1 = first_copy_train.cline

        # 2. combine T4 and R41, create Q4 at n4, R4' between n3 and n4
        pulse_highlight_box(VGroup(t4_down["shape"], r41), pulse_count=2, pulse_time=0.16)
        n4_w3_merge = merge_wire_tt_block(
            wire_y=main_y(3),
            hd_gate=None,
            target_gates=[t4_down["shape"], r41],
            ti_label="temp",
            ti_color=COLOR_T_US,
            show_label=False,
            attached_leg_lines=[bond34_down, copy_cline1]
        )
        # Split -> V4 (stays at wire 3), R4' (moves to wire 2)
        n4_v4_split = split_tensor_to_right_canon(
            tensor_info=n4_w3_merge,
            current_wire_y=main_y(3),
            neighbor_wire_y=main_y(2),
            right_canon_label="Q_{4}",
            remainder_label="R_{4}'",
            right_canon_color=COLOR_PURPLE,
            remainder_color=COLOR_T_US,
            attached_leg_lines=[bond34_down, copy_cline1],
            move_time=0.75
        )

        # 3. combine U3 and R31, and merge it with R4' from below and name it R3
        pulse_highlight_box(VGroup(u3_down["shape"], r31, n4_v4_split["remainder_tensor"][0]), pulse_count=2, pulse_time=0.16)
        n3_w2_merge = merge_wire_tt_block(
            wire_y=main_y(2),
            hd_gate=None,
            target_gates=[u3_down["shape"], r31],
            ti_label="temp",
            ti_color=COLOR_T_US,
            show_label=False,
            attached_leg_lines=[bond23_down] # cline1 is coming from n4_v4_split
        )
        n3_r3_bridge = merge_bridge_tensor_with_wire_tensor(
            right_canon_split_info=n4_v4_split,
            wire_tensor_info=n3_w2_merge,
            output_label="R_{3}",
            show_label=True,
            move_time=0.45,
            fuse_time=0.45,
            remove_leg_lines=[bond34_down],
            extra_leg_lines=[bond23_down, copy_cline1]
        )

        # 4. Split R3 -> V3 and R3' [moving to wire 1]
        n3_v3_split = split_tensor_to_right_canon(
            tensor_info=n3_r3_bridge,
            current_wire_y=main_y(2),
            neighbor_wire_y=main_y(1),
            right_canon_label="Q_{3}",
            remainder_label="R_{3}'",
            right_canon_color=COLOR_PURPLE,
            remainder_color=COLOR_T_US,
            attached_leg_lines=[bond23_down, copy_cline1],
            move_time=0.72
        )

        # 5. combine U2 and R21, and merge with R3' -> R2
        pulse_highlight_box(VGroup(u2_down["shape"], r21, n3_v3_split["remainder_tensor"][0]), pulse_count=2, pulse_time=0.16)
        n2_w1_merge = merge_wire_tt_block(
            wire_y=main_y(1),
            hd_gate=None,
            target_gates=[u2_down["shape"], r21],
            ti_label="temp",
            ti_color=COLOR_T_US,
            show_label=False,
            attached_leg_lines=[bond12_down]
        )
        n2_r2_bridge = merge_bridge_tensor_with_wire_tensor(
            right_canon_split_info=n3_v3_split,
            wire_tensor_info=n2_w1_merge,
            output_label="R_{2}",
            show_label=True,
            move_time=0.45,
            fuse_time=0.45,
            remove_leg_lines=[bond23_down],
            extra_leg_lines=[bond12_down, copy_cline1]
        )

        self.wait(1.0)

        # ================================================================
        # PHASE 6: ENTRY INTO COPY REGISTER n1'
        # ================================================================

        # --- Step 1: Expand control dot to wide R1' and shift n1' to midpoint ---
        pulse_highlight_box(VGroup(n2_r2_bridge["shape"], copy_dot1), pulse_count=2, pulse_time=0.16)
        r2_box = n2_r2_bridge["shape"][2]
        r2_x_left, r2_x_right = r2_box.get_left()[0], r2_box.get_right()[0]
        
        # Current y and target y
        y_n1 = main_y(0)
        y_n2 = main_y(1)
        midpoint_y = (y_n1 + y_n2) / 2
        old_copy_y0 = copy_y(0)
        
        # Store control line's original x position (must do BEFORE transforming copy_dot1)
        control_x = copy_dot1.get_center()[0]
        
        r1_prime_shape = make_wide_hatched_tensor(r2_x_left, r2_x_right, copy_y(0), COLOR_T_US)
        r1_prime_label = MathTex("R_{1}'", color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE).move_to(r1_prime_shape[2].get_center()).set_z_index(8)
        
        # Expansion
        self.play(
            ReplacementTransform(copy_dot1, r1_prime_shape),
            Write(r1_prime_label),
            run_time=0.6
        )

        # Move n1' wire and R1' to midpoint, keeping control line vertical
        copy_wire_1 = wires[1]
        n1p_shift = midpoint_y - old_copy_y0
        step1_anims = [
            copy_wire_1.animate.set_y(midpoint_y),
            r1_prime_shape.animate.set_y(midpoint_y),
            r1_prime_label.animate.set_y(midpoint_y),
            copy_cline1.animate.put_start_and_end_on(
                RIGHT * control_x + UP * midpoint_y,
                copy_cline1.get_end()
            ),
        ]
        self.play(*step1_anims, run_time=0.7)
        # Update logical position tracking
        wire_y_positions[1] = midpoint_y
        new_copy_y0 = midpoint_y

        # --- Step 2a: First, expand space to give room for QR ---
        y_cut = (new_copy_y0 + main_y(1)) / 2
        spacing_delta = 0.5  # Expand by this amount
        
        expand_anims = adjust_circuit_spacing(y_cut, spacing_delta, duration=0.5)
        self.play(*expand_anims)
        
        # After expansion, wires will have moved by delta/2
        n2_y_after_expansion = main_y(1) - spacing_delta / 2
        n1p_y_after_expansion = new_copy_y0 + spacing_delta / 2
        
        # --- Step 2b: Perform QR split using the established pattern ---
        # Attach label to shape for proper tracking
        n2_r2_bridge["shape"].label = n2_r2_bridge["label"]
        
        # Split R2 into Q2 (at n2) and R2' (at midpoint)
        n2_r2_split = split_tensor_to_right_canon(
            tensor_info=n2_r2_bridge,
            current_wire_y=n2_y_after_expansion,
            neighbor_wire_y=n1p_y_after_expansion,
            right_canon_label="Q_{2}",
            remainder_label="R_{2}'",
            attached_leg_lines=[bond12_down, copy_cline1],
            right_canon_color=COLOR_PURPLE,
            remainder_color=COLOR_T_US,
            move_time=0.75
        )
        
        # --- Step 3: Merge R2' with R1' to create T1' ---
        # Prepare wire tensor info for R1'
        r1_prime_info = {
            "shape": r1_prime_shape,
            "label": r1_prime_label,
        }
        
        # Attach label to shape for proper tracking
        r1_prime_shape.label = r1_prime_label
        
        # Merge R2' from above with R1' on the wire
        # Remove copy_cline1 (consumed in merge), but keep bond12_down (connects U1 to T1')
        n1_prime_bridge = merge_bridge_tensor_with_wire_tensor(
            right_canon_split_info=n2_r2_split,
            wire_tensor_info=r1_prime_info,
            output_label="temp",
            show_label=False,
            move_time=0.45,
            fuse_time=0.45,
            remove_leg_lines=[copy_cline1],
            extra_leg_lines=[]
        )
        
        # --- Step 4: Contract space back to original ---
        contract_anims = adjust_circuit_spacing(y_cut, -spacing_delta, duration=0.5)
        self.play(*contract_anims)
        
        # --- Step 5: Reshape to square and align U1 ---
        merged_tensor = n1_prime_bridge["shape"]
        x_left = merged_tensor[2].get_left()[0]
        x_right = merged_tensor[2].get_right()[0]
        pivot_x = (x_left + x_right) / 2
        y_merge = merged_tensor[2].get_center()[1]
        
        # Create square T1' as target
        square_t1p_shape = make_wide_hatched_tensor(
            pivot_x - 0.5 * GATE_WIDTH,
            pivot_x + 0.5 * GATE_WIDTH,
            y_merge,
            COLOR_T_US
        )
        t1_prime_label = MathTex("T_{1}'", color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
        t1_prime_label.move_to(square_t1p_shape[2].get_center()).set_z_index(8)
        
        # Get SVD bond and Q2 from split/merge results
        svd_bond = n1_prime_bridge["svd_bond"]
        q2_shape = n2_r2_split["right_canon_tensor"][0]
        q2_label = n2_r2_split["right_canon_tensor"][1]
        
        # Add updaters to bonds so they track during reshape
        svd_bond.add_updater(
            lambda mob: mob.put_start_and_end_on(
                q2_shape[2].get_center(),
                square_t1p_shape[2].get_top()
            )
        )
        bond12_down.add_updater(
            lambda mob: mob.put_start_and_end_on(
                np.array([u1_down["shape"][2].get_center()[0], u1_down["shape"][2].get_top()[1], 0.0]),
                np.array([u1_down["shape"][2].get_center()[0], square_t1p_shape[2].get_top()[1], 0.0]),
            )
        )
        
        # Transform merged tensor to square, align U1
        self.play(
            ReplacementTransform(merged_tensor, square_t1p_shape),
            Write(t1_prime_label),
            u1_down["shape"].animate.set_x(pivot_x),
            u1_down["label"].animate.set_x(pivot_x),
            run_time=0.6
        )
        
        # Clear updaters and fix final positions
        svd_bond.clear_updaters()
        bond12_down.clear_updaters()
        svd_bond.put_start_and_end_on(q2_shape[2].get_center(), square_t1p_shape[2].get_top())
        bond12_down.put_start_and_end_on(
            np.array([u1_down["shape"][2].get_center()[0], u1_down["shape"][2].get_top()[1], 0.0]),
            np.array([u1_down["shape"][2].get_center()[0], square_t1p_shape[2].get_top()[1], 0.0]),
        )
        self.wait(0.6)
        
        # Store result for downstream use
        t1_prime_shape = square_t1p_shape
        n1_prime_merge = {
            "tensor": VGroup(t1_prime_shape, t1_prime_label),
            "shape": t1_prime_shape,
            "label": t1_prime_label,
            "svd_bond": svd_bond,
            "right_canon_tensor": VGroup(q2_shape, q2_label),
            "legs": [bond12_down],
        }

        # --- Step 5: Downward recursive sweep producing U' registers ---
        # Map Q2, Q3, Q4 for the recursion
        phase6_canon_map = {
            2: {"shape": n2_r2_split["right_canon_tensor"][0], "label": n2_r2_split["right_canon_tensor"][1], "svd_bond": n2_r2_split.get("svd_bond")},
            3: {"shape": n3_v3_split["right_canon_tensor"][0], "label": n3_v3_split["right_canon_tensor"][1], "svd_bond": n3_v3_split.get("svd_bond")},
            4: {"shape": n4_v4_split["right_canon_tensor"][0], "label": n4_v4_split["right_canon_tensor"][1], "svd_bond": n4_v4_split.get("svd_bond")},
        }

        # Perform the recursion downward through registers [n1', n2, n3, n4]
        final_t_p6, last_u_p6, last_bond_p6 = compress_oc_recursive(
            train_idx=1,
            t_current=n1_prime_merge,
            canon_map=phase6_canon_map,
            direction="down",
            previous_canon=u1_down,
            blue_bond_to_current_t=bond12_down,
            show_pulse=False,
            y_override_list=[new_copy_y0, main_y(1), main_y(2), main_y(3)],
            label_suffix="'",
            prime_canon_indices={1},
        )
        self.wait(0.6)

        # ================================================================
        # PHASE 7: ENTRY INTO COPY REGISTER n2'
        # ================================================================
        # After Phase 6, we have:
        # - U1' at n1', U2 at n2, U3 at n3, T4 at n4 (final_t_p6)
        # - Second copy train: copy_items[1] has dot at n2', gates R32 at n3, R42 at n4
        
        # Extract components from Phase 6 result
        u3_p6 = last_u_p6
        bond34_p6 = last_bond_p6
        u2_p6 = u3_p6["upstream_canon"]
        bond23_p6 = u3_p6["upstream_bond"]
        u1_prime_p6 = u2_p6["upstream_canon"]
        bond12_p6 = u2_p6["upstream_bond"]
        u1_main_p6 = u1_prime_p6["upstream_canon"] if u1_prime_p6 is not None else None
        bond_u1_to_u1p_p6 = u1_prime_p6["upstream_bond"] if u1_prime_p6 is not None else None
        t4_p6 = final_t_p6

        # Access second copy train
        second_copy_train = copy_items[1]
        r32 = second_copy_train.gates[0]  # R32 at n3
        r42 = second_copy_train.gates[1]  # R42 at n4
        copy_dot2 = second_copy_train.dot
        copy_cline2 = second_copy_train.cline

        # Attach labels to shapes for merge handling
        t4_p6["shape"].label = t4_p6["label"]
        u3_p6["shape"].label = u3_p6["label"]

        # --- Step 1: Merge T4 with R42 -> wide tensor, then QR split ---
        pulse_highlight_box(VGroup(t4_p6["shape"], r42), pulse_count=2, pulse_time=0.16)
        n4_w3_merge_p7 = merge_wire_tt_block(
            wire_y=main_y(3),
            hd_gate=None,
            target_gates=[t4_p6["shape"], r42],
            ti_label="temp",
            ti_color=COLOR_T_US,
            show_label=False,
            attached_leg_lines=[bond34_p6, copy_cline2]
        )
        n4_v4_split_p7 = split_tensor_to_right_canon(
            tensor_info=n4_w3_merge_p7,
            current_wire_y=main_y(3),
            neighbor_wire_y=main_y(2),
            right_canon_label="Q_{4}",
            remainder_label="R_{4}'",
            right_canon_color=COLOR_PURPLE,
            remainder_color=COLOR_T_US,
            attached_leg_lines=[bond34_p6, copy_cline2],
            move_time=0.75
        )

        # --- Step 2: Merge U3 with R32, then merge R4' from above -> R3 ---
        pulse_highlight_box(VGroup(u3_p6["shape"], r32, n4_v4_split_p7["remainder_tensor"][0]), pulse_count=2, pulse_time=0.16)
        n3_w2_merge_p7 = merge_wire_tt_block(
            wire_y=main_y(2),
            hd_gate=None,
            target_gates=[u3_p6["shape"], r32],
            ti_label="temp",
            ti_color=COLOR_T_US,
            show_label=False,
            attached_leg_lines=[bond23_p6]
        )
        n3_r3_bridge_p7 = merge_bridge_tensor_with_wire_tensor(
            right_canon_split_info=n4_v4_split_p7,
            wire_tensor_info=n3_w2_merge_p7,
            output_label="R_{3}",
            show_label=True,
            move_time=0.45,
            fuse_time=0.45,
            remove_leg_lines=[bond34_p6],
            extra_leg_lines=[bond23_p6, copy_cline2]
        )

        # --- Step 3: Expand control dot to wide R2' and shift n2' to midpoint ---
        pulse_highlight_box(VGroup(n3_r3_bridge_p7["shape"], copy_dot2), pulse_count=2, pulse_time=0.16)
        y_n2 = main_y(1)
        y_n3 = main_y(2)
        midpoint_y_p7 = (y_n2 + y_n3) / 2
        old_copy_y1 = copy_y(1)
        
        control_x_p7 = copy_dot2.get_center()[0]
        r3_box_p7 = n3_r3_bridge_p7["shape"][2]
        r3_x_left_p7, r3_x_right_p7 = r3_box_p7.get_left()[0], r3_box_p7.get_right()[0]
        
        r2_prime_shape_p7 = make_wide_hatched_tensor(r3_x_left_p7, r3_x_right_p7, copy_y(1), COLOR_T_US)
        r2_prime_label_p7 = MathTex("R_{2}'", color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE).move_to(r2_prime_shape_p7[2].get_center()).set_z_index(8)
        
        self.play(
            ReplacementTransform(copy_dot2, r2_prime_shape_p7),
            Write(r2_prime_label_p7),
            run_time=0.6
        )

        # Move n2' wire and R2' to midpoint
        copy_wire_2 = wires[3]  # n2' is wire index 3
        n2p_shift = midpoint_y_p7 - old_copy_y1
        step1_anims_p7 = [
            copy_wire_2.animate.set_y(midpoint_y_p7),
            r2_prime_shape_p7.animate.set_y(midpoint_y_p7),
            r2_prime_label_p7.animate.set_y(midpoint_y_p7),
            copy_cline2.animate.put_start_and_end_on(
                RIGHT * control_x_p7 + UP * midpoint_y_p7,
                copy_cline2.get_end()
            ),
        ]
        self.play(*step1_anims_p7, run_time=0.7)
        wire_y_positions[3] = midpoint_y_p7
        new_copy_y1 = midpoint_y_p7

        # --- Step 4: Expand space, QR split, merge, contract ---
        y_cut_p7 = (new_copy_y1 + main_y(2)) / 2
        spacing_delta_p7 = 0.5
        
        expand_anims_p7 = adjust_circuit_spacing(y_cut_p7, spacing_delta_p7, duration=0.5)
        self.play(*expand_anims_p7)
        
        n3_y_after_expansion_p7 = main_y(2) - spacing_delta_p7 / 2
        n2p_y_after_expansion_p7 = new_copy_y1 + spacing_delta_p7 / 2
        
        # QR split R3 -> Q3, R3'
        n3_r3_bridge_p7["shape"].label = n3_r3_bridge_p7["label"]
        n3_r3_split_p7 = split_tensor_to_right_canon(
            tensor_info=n3_r3_bridge_p7,
            current_wire_y=n3_y_after_expansion_p7,
            neighbor_wire_y=n2p_y_after_expansion_p7,
            right_canon_label="Q_{3}",
            remainder_label="R_{3}'",
            attached_leg_lines=[bond23_p6, copy_cline2],
            right_canon_color=COLOR_PURPLE,
            remainder_color=COLOR_T_US,
            move_time=0.75
        )
        
        # Merge R3' with R2'
        r2_prime_info_p7 = {
            "shape": r2_prime_shape_p7,
            "label": r2_prime_label_p7,
        }
        r2_prime_shape_p7.label = r2_prime_label_p7
        
        n2_prime_bridge_p7 = merge_bridge_tensor_with_wire_tensor(
            right_canon_split_info=n3_r3_split_p7,
            wire_tensor_info=r2_prime_info_p7,
            output_label="temp",
            show_label=False,
            move_time=0.45,
            fuse_time=0.45,
            remove_leg_lines=[copy_cline2],
            extra_leg_lines=[]
        )
        
        # Contract space
        contract_anims_p7 = adjust_circuit_spacing(y_cut_p7, -spacing_delta_p7, duration=0.5)
        self.play(*contract_anims_p7)
        
        # Reshape to square T2' and align U2
        merged_tensor_p7 = n2_prime_bridge_p7["shape"]
        x_left_p7 = merged_tensor_p7[2].get_left()[0]
        x_right_p7 = merged_tensor_p7[2].get_right()[0]
        pivot_x_p7 = (x_left_p7 + x_right_p7) / 2
        y_merge_p7 = merged_tensor_p7[2].get_center()[1]
        
        square_t2p_shape = make_wide_hatched_tensor(
            pivot_x_p7 - 0.5 * GATE_WIDTH,
            pivot_x_p7 + 0.5 * GATE_WIDTH,
            y_merge_p7,
            COLOR_T_US
        )
        t2_prime_label = MathTex("T_{2}'", color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
        t2_prime_label.move_to(square_t2p_shape[2].get_center()).set_z_index(8)
        
        svd_bond_p7 = n2_prime_bridge_p7["svd_bond"]
        q3_shape_p7 = n3_r3_split_p7["right_canon_tensor"][0]
        q3_label_p7 = n3_r3_split_p7["right_canon_tensor"][1]
        
        svd_bond_p7.add_updater(
            lambda mob: mob.put_start_and_end_on(
                q3_shape_p7[2].get_center(),
                square_t2p_shape[2].get_top()
            )
        )
        bond23_p6.add_updater(
            lambda mob: mob.put_start_and_end_on(
                np.array([u2_p6["shape"][2].get_center()[0], u2_p6["shape"][2].get_top()[1], 0.0]),
                np.array([u2_p6["shape"][2].get_center()[0], square_t2p_shape[2].get_top()[1], 0.0]),
            )
        )
        
        u2_old_x_p7 = u2_p6["shape"][2].get_center()[0]
        pivot_shift_x_p7 = pivot_x_p7 - u2_old_x_p7
        pivot_anims_p7 = [
            ReplacementTransform(merged_tensor_p7, square_t2p_shape),
            Write(t2_prime_label),
            u2_p6["shape"].animate.set_x(pivot_x_p7),
            u2_p6["label"].animate.set_x(pivot_x_p7),
            u1_prime_p6["shape"].animate.set_x(pivot_x_p7),
            u1_prime_p6["label"].animate.set_x(pivot_x_p7),
            bond12_p6.animate.shift(RIGHT * pivot_shift_x_p7),
        ]
        if u1_main_p6 is not None:
            pivot_anims_p7.extend(
                [
                    u1_main_p6["shape"].animate.set_x(pivot_x_p7),
                    u1_main_p6["label"].animate.set_x(pivot_x_p7),
                ]
            )
        if bond_u1_to_u1p_p6 is not None:
            pivot_anims_p7.append(bond_u1_to_u1p_p6.animate.shift(RIGHT * pivot_shift_x_p7))
        self.play(*pivot_anims_p7, run_time=0.6)
        
        svd_bond_p7.clear_updaters()
        bond23_p6.clear_updaters()
        svd_bond_p7.put_start_and_end_on(q3_shape_p7[2].get_center(), square_t2p_shape[2].get_top())
        bond23_p6.put_start_and_end_on(
            np.array([u2_p6["shape"][2].get_center()[0], u2_p6["shape"][2].get_top()[1], 0.0]),
            np.array([u2_p6["shape"][2].get_center()[0], square_t2p_shape[2].get_top()[1], 0.0]),
        )
        self.wait(0.6)
        
        t2_prime_shape = square_t2p_shape
        n2_prime_merge_p7 = {
            "tensor": VGroup(t2_prime_shape, t2_prime_label),
            "shape": t2_prime_shape,
            "label": t2_prime_label,
            "svd_bond": svd_bond_p7,
            "right_canon_tensor": VGroup(q3_shape_p7, q3_label_p7),
            "legs": [bond23_p6],
        }

        # Downward recursive sweep
        phase7_canon_map = {
            3: {"shape": n3_r3_split_p7["right_canon_tensor"][0], "label": n3_r3_split_p7["right_canon_tensor"][1], "svd_bond": n3_r3_split_p7.get("svd_bond")},
            4: {"shape": n4_v4_split_p7["right_canon_tensor"][0], "label": n4_v4_split_p7["right_canon_tensor"][1], "svd_bond": n4_v4_split_p7.get("svd_bond")},
        }

        final_t_p7, last_u_p7, last_bond_p7 = compress_oc_recursive(
            train_idx=2,
            t_current=n2_prime_merge_p7,
            canon_map=phase7_canon_map,
            direction="down",
            previous_canon=u2_p6,
            blue_bond_to_current_t=bond23_p6,
            show_pulse=False,
            y_override_list=[new_copy_y0, new_copy_y1, main_y(2), main_y(3)],
            label_suffix="'",
            prime_canon_indices={2},
        )
        self.wait(0.6)

        # ================================================================
        # PHASE 8: ENTRY INTO COPY REGISTER n3'
        # ================================================================
        # After Phase 7, we have:
        # - U2' at n2', U3 at n3, T4 at n4 (final_t_p7)
        # - Third copy train: copy_items[2] has dot at n3', gate R43 at n4
        
        u3_p7 = last_u_p7
        bond34_p7 = last_bond_p7
        u2_prime_p7 = u3_p7["upstream_canon"]
        bond23_p7 = u3_p7["upstream_bond"]
        u2_main_p7 = u2_prime_p7["upstream_canon"] if u2_prime_p7 is not None else None
        bond_u2_to_u2p_p7 = u2_prime_p7["upstream_bond"] if u2_prime_p7 is not None else None
        u1_prime_p7 = u2_main_p7["upstream_canon"] if u2_main_p7 is not None else None
        bond_u1p_to_u2_p7 = u2_main_p7["upstream_bond"] if u2_main_p7 is not None else None
        u1_main_p7 = u1_prime_p7["upstream_canon"] if u1_prime_p7 is not None else None
        bond_u1_to_u1p_p7 = u1_prime_p7["upstream_bond"] if u1_prime_p7 is not None else None
        t4_p7 = final_t_p7

        # Access third copy train
        third_copy_train = copy_items[2]
        r43 = third_copy_train.gates[0]  # R43 at n4
        copy_dot3 = third_copy_train.dot
        copy_cline3 = third_copy_train.cline

        # Attach labels
        t4_p7["shape"].label = t4_p7["label"]
        u3_p7["shape"].label = u3_p7["label"]

        # --- Step 1: Merge T4 with R43 -> R4 on n4 ---
        pulse_highlight_box(VGroup(t4_p7["shape"], r43), pulse_count=2, pulse_time=0.16)
        n4_r4_bridge_p8 = merge_wire_tt_block(
            wire_y=main_y(3),
            hd_gate=None,
            target_gates=[t4_p7["shape"], r43],
            ti_label="R_{4}",
            ti_color=COLOR_T_US,
            show_label=True,
            attached_leg_lines=[bond34_p7, copy_cline3]
        )

        # --- Step 2: Expand control dot to wide R3' and shift n3' to midpoint ---
        pulse_highlight_box(VGroup(n4_r4_bridge_p8["shape"], copy_dot3), pulse_count=2, pulse_time=0.16)
        y_n3_p8 = main_y(2)
        y_n4 = main_y(3)
        midpoint_y_p8 = (y_n3_p8 + y_n4) / 2
        old_copy_y2 = copy_y(2)
        
        control_x_p8 = copy_dot3.get_center()[0]
        r4_box_p8 = n4_r4_bridge_p8["shape"][2]
        r4_x_left_p8, r4_x_right_p8 = r4_box_p8.get_left()[0], r4_box_p8.get_right()[0]
        
        r3_prime_shape_p8 = make_wide_hatched_tensor(r4_x_left_p8, r4_x_right_p8, copy_y(2), COLOR_T_US)
        r3_prime_label_p8 = MathTex("R_{3}'", color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE).move_to(r3_prime_shape_p8[2].get_center()).set_z_index(8)
        
        self.play(
            ReplacementTransform(copy_dot3, r3_prime_shape_p8),
            Write(r3_prime_label_p8),
            run_time=0.6
        )

        # Move n3' wire and R3' to midpoint
        copy_wire_3 = wires[5]  # n3' is wire index 5
        n3p_shift = midpoint_y_p8 - old_copy_y2
        step1_anims_p8 = [
            copy_wire_3.animate.set_y(midpoint_y_p8),
            r3_prime_shape_p8.animate.set_y(midpoint_y_p8),
            r3_prime_label_p8.animate.set_y(midpoint_y_p8),
            copy_cline3.animate.put_start_and_end_on(
                RIGHT * control_x_p8 + UP * midpoint_y_p8,
                copy_cline3.get_end()
            ),
        ]
        self.play(*step1_anims_p8, run_time=0.7)
        wire_y_positions[5] = midpoint_y_p8
        new_copy_y2 = midpoint_y_p8

        # --- Step 3: Expand space, QR split R4, merge into R3', contract ---
        y_cut_p8 = (new_copy_y2 + main_y(3)) / 2
        spacing_delta_p8 = 0.5
        
        expand_anims_p8 = adjust_circuit_spacing(y_cut_p8, spacing_delta_p8, duration=0.5)
        self.play(*expand_anims_p8)
        
        n4_y_after_expansion_p8 = main_y(3) - spacing_delta_p8 / 2
        n3p_y_after_expansion_p8 = new_copy_y2 + spacing_delta_p8 / 2
        
        # QR split R4 -> Q4, R4' after expansion
        n4_r4_bridge_p8["shape"].label = n4_r4_bridge_p8["label"]
        n4_v4_split_p8 = split_tensor_to_right_canon(
            tensor_info=n4_r4_bridge_p8,
            current_wire_y=n4_y_after_expansion_p8,
            neighbor_wire_y=n3p_y_after_expansion_p8,
            right_canon_label="Q_{4}",
            remainder_label="R_{4}'",
            right_canon_color=COLOR_PURPLE,
            remainder_color=COLOR_T_US,
            attached_leg_lines=[bond34_p7, copy_cline3],
            move_time=0.75
        )
        
        # Merge R4' with R3'
        r3_prime_info_p8 = {
            "shape": r3_prime_shape_p8,
            "label": r3_prime_label_p8,
        }
        r3_prime_shape_p8.label = r3_prime_label_p8
        
        n3_prime_bridge_p8 = merge_bridge_tensor_with_wire_tensor(
            right_canon_split_info=n4_v4_split_p8,
            wire_tensor_info=r3_prime_info_p8,
            output_label="temp",
            show_label=False,
            move_time=0.45,
            fuse_time=0.45,
            remove_leg_lines=[copy_cline3],
            extra_leg_lines=[]
        )
        
        # Contract space
        contract_anims_p8 = adjust_circuit_spacing(y_cut_p8, -spacing_delta_p8, duration=0.5)
        self.play(*contract_anims_p8)
        
        # Reshape to square T3' and align U3
        merged_tensor_p8 = n3_prime_bridge_p8["shape"]
        x_left_p8 = merged_tensor_p8[2].get_left()[0]
        x_right_p8 = merged_tensor_p8[2].get_right()[0]
        pivot_x_p8 = (x_left_p8 + x_right_p8) / 2
        y_merge_p8 = merged_tensor_p8[2].get_center()[1]
        
        square_t3p_shape = make_wide_hatched_tensor(
            pivot_x_p8 - 0.5 * GATE_WIDTH,
            pivot_x_p8 + 0.5 * GATE_WIDTH,
            y_merge_p8,
            COLOR_T_US
        )
        t3_prime_label = MathTex("T_{3}'", color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(GATE_LABEL_SCALE)
        t3_prime_label.move_to(square_t3p_shape[2].get_center()).set_z_index(8)
        
        svd_bond_p8 = n3_prime_bridge_p8["svd_bond"]
        q4_shape_p8 = n4_v4_split_p8["right_canon_tensor"][0]
        q4_label_p8 = n4_v4_split_p8["right_canon_tensor"][1]
        
        svd_bond_p8.add_updater(
            lambda mob: mob.put_start_and_end_on(
                q4_shape_p8[2].get_center(),
                square_t3p_shape[2].get_top()
            )
        )
        bond34_p7.add_updater(
            lambda mob: mob.put_start_and_end_on(
                np.array([u3_p7["shape"][2].get_center()[0], u3_p7["shape"][2].get_top()[1], 0.0]),
                np.array([u3_p7["shape"][2].get_center()[0], square_t3p_shape[2].get_top()[1], 0.0]),
            )
        )
        
        u3_old_x_p8 = u3_p7["shape"][2].get_center()[0]
        pivot_shift_x_p8 = pivot_x_p8 - u3_old_x_p8
        pivot_anims_p8 = [
            ReplacementTransform(merged_tensor_p8, square_t3p_shape),
            Write(t3_prime_label),
            u3_p7["shape"].animate.set_x(pivot_x_p8),
            u3_p7["label"].animate.set_x(pivot_x_p8),
            u2_prime_p7["shape"].animate.set_x(pivot_x_p8),
            u2_prime_p7["label"].animate.set_x(pivot_x_p8),
            bond23_p7.animate.shift(RIGHT * pivot_shift_x_p8),
        ]
        if u2_main_p7 is not None:
            pivot_anims_p8.extend(
                [
                    u2_main_p7["shape"].animate.set_x(pivot_x_p8),
                    u2_main_p7["label"].animate.set_x(pivot_x_p8),
                ]
            )
        if bond_u2_to_u2p_p7 is not None:
            pivot_anims_p8.append(bond_u2_to_u2p_p7.animate.shift(RIGHT * pivot_shift_x_p8))
        if u1_prime_p7 is not None:
            pivot_anims_p8.extend(
                [
                    u1_prime_p7["shape"].animate.set_x(pivot_x_p8),
                    u1_prime_p7["label"].animate.set_x(pivot_x_p8),
                ]
            )
        if bond_u1p_to_u2_p7 is not None:
            pivot_anims_p8.append(bond_u1p_to_u2_p7.animate.shift(RIGHT * pivot_shift_x_p8))
        if u1_main_p7 is not None:
            pivot_anims_p8.extend(
                [
                    u1_main_p7["shape"].animate.set_x(pivot_x_p8),
                    u1_main_p7["label"].animate.set_x(pivot_x_p8),
                ]
            )
        if bond_u1_to_u1p_p7 is not None:
            pivot_anims_p8.append(bond_u1_to_u1p_p7.animate.shift(RIGHT * pivot_shift_x_p8))
        self.play(*pivot_anims_p8, run_time=0.6)
        
        svd_bond_p8.clear_updaters()
        bond34_p7.clear_updaters()
        svd_bond_p8.put_start_and_end_on(q4_shape_p8[2].get_center(), square_t3p_shape[2].get_top())
        bond34_p7.put_start_and_end_on(
            np.array([u3_p7["shape"][2].get_center()[0], u3_p7["shape"][2].get_top()[1], 0.0]),
            np.array([u3_p7["shape"][2].get_center()[0], square_t3p_shape[2].get_top()[1], 0.0]),
        )
        self.wait(0.6)
        
        t3_prime_shape = square_t3p_shape
        n3_prime_merge_p8 = {
            "tensor": VGroup(t3_prime_shape, t3_prime_label),
            "shape": t3_prime_shape,
            "label": t3_prime_label,
            "svd_bond": svd_bond_p8,
            "right_canon_tensor": VGroup(q4_shape_p8, q4_label_p8),
            "legs": [bond34_p7],
        }

        # Final downward sweep to n4
        phase8_canon_map = {
            4: {"shape": n4_v4_split_p8["right_canon_tensor"][0], "label": n4_v4_split_p8["right_canon_tensor"][1], "svd_bond": n4_v4_split_p8.get("svd_bond")},
        }

        final_t_p8, last_u_p8, last_bond_p8 = compress_oc_recursive(
            train_idx=3,
            t_current=n3_prime_merge_p8,
            canon_map=phase8_canon_map,
            direction="down",
            previous_canon=u3_p7,
            blue_bond_to_current_t=bond34_p7,
            show_pulse=False,
            y_override_list=[new_copy_y0, new_copy_y1, new_copy_y2, main_y(3)],
            label_suffix="'",
            prime_canon_indices={3},
        )
        self.wait(0.6)

        # PHASE 9: IDENTITY COMPLETION AT n4' + MPO FORMATION
        # ================================================================

        # Move n4' so all main-copy spacings are uniform after midpoint moves.
        old_copy_y3 = copy_y(3)
        target_copy_y3 = main_y(3) - 0.5 * (main_y(2) - main_y(3))
        n4p_shift = target_copy_y3 - old_copy_y3
        move_n4p_anims = [wires[7].animate.set_y(target_copy_y3)]
        self.play(*move_n4p_anims, run_time=0.6)
        wire_y_positions[7] = target_copy_y3

        # Create damping identity at n4'
        damp_pivot_x = final_t_p8["shape"][2].get_center()[0]
        damping_i_n4p = QuantumGate("I", damp_pivot_x, target_copy_y3, WHITE)
        damping_i_n4p.box.set_z_index(5)
        damping_i_n4p.label.set_z_index(8)
        self.play(
            Create(damping_i_n4p.box),
            Write(damping_i_n4p.label),
            run_time=0.7,
        )
        self.add(damping_i_n4p)
        
        # Attach damping identity bond from n4 to n4'
        damping_i_bond = Line(
            final_t_p8["shape"][2].get_bottom(),
            damping_i_n4p.box.get_top(),
            color=COLOR_LIGHT_BLUE,
            stroke_width=BOND_LIGHT,
        ).set_z_index(0)
        self.play(Create(damping_i_bond), run_time=0.45)

        self.wait(1.0)

        # Gather all MPO components to shift together (Don't fuse, keep labels/bonds)
        mpo_mobs = []
        
        # Follow the chain from last_u_p8 backwards
        cursor = last_u_p8
        while cursor is not None:
            mpo_mobs.append(cursor["shape"])
            mpo_mobs.append(cursor["label"])
            if "upstream_bond" in cursor and cursor["upstream_bond"] is not None:
                mpo_mobs.append(cursor["upstream_bond"])
            cursor = cursor.get("upstream_canon")
        
        # Add final_t_p8, damping identity, and their connecting bonds
        mpo_mobs.append(final_t_p8["shape"])
        mpo_mobs.append(final_t_p8["label"])
        mpo_mobs.append(last_bond_p8)
        mpo_mobs.append(damping_i_n4p.box)
        mpo_mobs.append(damping_i_n4p.label)
        mpo_mobs.append(damping_i_bond)
        
        mpo_group = VGroup(*[mob for mob in mpo_mobs if mob is not None])
        
        # Center final MPO and contract wires to short MPO legs
        recenter_dx = CIRCUIT_CENTER_X - damp_pivot_x
        wire_leg_half = 0.6  # Slightly longer legs
        
        final_layout_anims = [mpo_group.animate.shift(RIGHT * recenter_dx)]
        
        # Shift residual vertical lines if any (like copy clines that weren't consumed?)
        protected_ids = {id(mob) for mob in mpo_group.get_family()}
        wire_family_ids = {id(m) for w in wires for m in w.get_family()}
        protected_ids.update(wire_family_ids)
        
        residual_mobs = [
            mob for mob in self.mobjects 
            if id(mob) not in protected_ids and isinstance(mob, (MathTex, Line))
        ]
        if abs(recenter_dx) > 1e-8:
            final_layout_anims.extend([mob.animate.shift(RIGHT * recenter_dx) for mob in residual_mobs])

        # Contract wires
        for w in wires:
            y = w.line.get_center()[1]
            new_start = np.array([CIRCUIT_CENTER_X - wire_leg_half, y, 0.0])
            new_end = np.array([CIRCUIT_CENTER_X + wire_leg_half, y, 0.0])
            final_layout_anims.extend([
                w.line.animate.put_start_and_end_on(new_start, new_end),
                w.start_circ.animate.move_to(new_start),
                w.end_circ.animate.move_to(new_end),
                w.label.animate.move_to(new_start + LEFT * 0.4),
            ])
            
        self.play(*final_layout_anims, run_time=1.2)
        self.wait(1.5)
        
