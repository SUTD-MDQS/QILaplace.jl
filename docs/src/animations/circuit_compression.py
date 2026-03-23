# If you are here looking for the source codes for my beautiful animation, a fair warning!
# This script is mostly vibe-coded using LLMs and they might not be the elegant script you are searching for. I only wrote codes upto line 254 in this one. The rest is stitched up by ChatGPT. Good luck. 
from manim import *
import shutil
import os
from pathlib import Path

config.background_color = "#222831"

class CircuitCompression(Scene):
    def construct(self):
        Text.set_default(color="#DDDDDD")
        MathTex.set_default(color="#DDDDDD")
        Tex.set_default(color="#DDDDDD")

        # ================================================================
        # PSEUDOCODE OVERVIEW
        # 1) Build visual primitives for a 4-qubit QFT-like circuit.
        # 2) Reveal the raw circuit (H gates + controlled phase trains).
        # 3) ZIP-UP sweep (bottom -> top): merge local gates into T_i and
        #    split into U_iS_i + V_i while updating wire spacing and bonds.
        # 4) PUSH-DOWN sweep (top -> bottom): recursively merge/split to
        #    orthogonalize and form the MPO chain.
        # 5) Run custom right-edge choreography (n4/n3 merge path), then
        #    continue recursion to complete compression.
        # 6) Recenter final tensors/wires and copy rendered media to assets.
        # ================================================================

        # -------------------------- Theme and style --------------------------
        # Theme Colors
        COLOR_PURPLE = "#8E47B7"      # V
        COLOR_GREEN = "#5F9EA0"       # T (Orthogonal Center - Dull Blue/CadetBlue)
        COLOR_LIGHT_BLUE = "#297DB5"  # Compressed bonds & S
        COLOR_ORANGE_RED = "#B14628"  # Heavy uncompressed bonds
        COLOR_BLACK = BLACK
        COLOR_T_US = "#87CEEB"        # Sky Blue for all T and U*S
        
        COLOR_H = "#3CB371"       # Greenish shade for H gate border
        COLOR_BROWN_1 = "#D98E43" # Darker brown (pi/2)
        COLOR_BROWN_2 = "#A47243" # Dark brown (pi/4)
        COLOR_BROWN_3 = "#6A4329" # Medium brown (pi/8)
        COLOR_YELLOW = "#EAD54B"  # U 
        
        # Configuration for gate dimensions and wire weight
        GATE_WIDTH = 1.0
        WIRE_THICKNESS = 6
        BOND_HEAVY = 18       # Width for uncompressed (orange-red) bonds
        BOND_LIGHT = 8        # Width for compressed (light-blue) bonds
        
        # Class representing a single qubit wire (ni)
        class QuantumWire(VGroup):
            def __init__(self, y_pos, length, start_x, idx):
                super().__init__()
                self.y_pos = y_pos
                self.idx = idx
                
                # Main horizontal line representing the qubit
                self.line = Line(
                    start=RIGHT * start_x + UP * y_pos,
                    end=RIGHT * (start_x + length) + UP * y_pos,
                    stroke_width=WIRE_THICKNESS,
                    color=WHITE
                )
                
                # Rounded endpoints for aesthetic appeal
                circ_rad = (WIRE_THICKNESS * 1.25) / 100 
                self.start_circ = Circle(radius=circ_rad, color=WHITE, fill_opacity=1).move_to(self.line.get_start())
                self.end_circ = Circle(radius=circ_rad, color=WHITE, fill_opacity=1).move_to(self.line.get_end())
                
                # Qubit index label (e.g., n1)
                self.label = MathTex(f"n_{idx}").next_to(self.start_circ, LEFT, buff=0.2)
                
                self.add(self.line, self.start_circ, self.end_circ, self.label)

        # Helper to create single gates (H)
        class QuantumGate(VGroup):
            def __init__(self, g_type, x_pos, target_y, color):
                super().__init__()
                self.x_pos = x_pos
                self.target_y = target_y
                
                side = GATE_WIDTH
                self.box = RoundedRectangle(corner_radius=0.1, width=side, height=side, 
                                            color=COLOR_BLACK, stroke_width=WIRE_THICKNESS * side, 
                                            fill_color=color, fill_opacity=1)

                self.box.move_to(RIGHT * x_pos + UP * target_y)
                self.label = MathTex(g_type, color=BLACK, stroke_color=BLACK, stroke_width=1.5).move_to(self.box.get_center())
                
                self.add(self.box, self.label)

        # Helper to create controlled phase trains
        # Class representing a group of phase gates controlled by a single wire
        class PhaseTrain(VGroup):
            def __init__(self, x_pos, control_y, targets_info):
                super().__init__()
                self.x_pos = x_pos
                self.control_y = control_y
                self.targets_info = targets_info
                self.gates = []
                
                # Control dot on the source wire
                dot_radius = (WIRE_THICKNESS * 2.5) / 100
                self.dot = Dot(point=RIGHT * x_pos + UP * control_y, radius=dot_radius, color=WHITE)
                
                # Vertical line connecting the control dot down to target gates
                min_y = min([t[0] for t in targets_info])
                self.cline = Line(self.dot.get_center(), RIGHT * x_pos + UP * min_y, stroke_width=WIRE_THICKNESS, color=WHITE).set_z_index(0)
                
                self.add(self.dot, self.cline)
                
                # Individual phase gates on target wires
                side = GATE_WIDTH
                for (y_pos, angle_str, color) in targets_info:
                    box = RoundedRectangle(corner_radius=0.1, width=side, height=side, 
                                           color=COLOR_BLACK, stroke_width=WIRE_THICKNESS * side, 
                                           fill_color=color, fill_opacity=1).set_z_index(1)
                    box.move_to(RIGHT * x_pos + UP * y_pos)
                    label_str = f"P(\\frac{{\\pi}}{{{angle_str}}})"
                    label = MathTex(label_str, color=BLACK, stroke_color=BLACK, stroke_width=1).move_to(box.get_center()).scale(0.65).set_z_index(2)
                    gate_group = VGroup(box, label)
                    self.gates.append(gate_group)
                    self.add(gate_group)

        class TensorNode(VGroup):
            def __init__(self, label_str, x_pos, y_pos, t_type="U"):
                super().__init__()
                color = COLOR_YELLOW if t_type == "U" else (COLOR_PURPLE if t_type == "V" else COLOR_GREEN)
                self.box = RoundedRectangle(corner_radius=0.1, width=1.0, height=1.0, 
                                            color=COLOR_BLACK, stroke_width=WIRE_THICKNESS * 1.0,
                                            fill_color=color, fill_opacity=1).move_to(RIGHT * x_pos + UP * y_pos)
                self.label = MathTex(label_str, color=BLACK, stroke_color=BLACK, stroke_width=1).move_to(self.box.get_center())
                self.add(self.box, self.label)
                self.set_z_index(5)

        # -------------------- Layout controls (variable glossary) --------------------
        # qubits: total number of wires in the animation.
        # MIN_GATE_GAP: minimum spacing between neighboring gate centers.
        # GATE_CENTER_SPACING: fixed pitch of timeline slots.
        # H_AND_TRAIN_COUNT: number of alternating [H, train, ..., H] slots.
        # WIRE_EDGE_MARGIN: left/right padding from wire endpoint to nearest gate slot.
        # CIRCUIT_CENTER_X: horizontal anchor used for intermediate and final recentering.
        # -----------------------------------------------------------------------------
        # Draw wires of the circuit
        qubits = 4
        MIN_GATE_GAP = 0.3
        GATE_CENTER_SPACING = GATE_WIDTH + MIN_GATE_GAP
        H_AND_TRAIN_COUNT = (2 * qubits) - 1
        WIRE_EDGE_MARGIN = MIN_GATE_GAP
        CIRCUIT_CENTER_X = 0.0

        circuit_wire_length = (
            (2 * WIRE_EDGE_MARGIN)
            + GATE_WIDTH
            + ((H_AND_TRAIN_COUNT - 1) * GATE_CENTER_SPACING)
        )
        wire_start_x = CIRCUIT_CENTER_X - (circuit_wire_length / 2)

        first_element_x = wire_start_x + WIRE_EDGE_MARGIN + (GATE_WIDTH / 2)
        element_x_positions = [
            first_element_x + i * GATE_CENTER_SPACING
            for i in range(H_AND_TRAIN_COUNT)
        ]

        # Initialization: Create the 4 qubit wires sequentially from top to bottom
        wire_spacing = 1.5
        wires = []
        for i in range(qubits):
            # n_1 is top-most wire, n_4 is bottom-most
            w = QuantumWire(
                y_pos=2 - i * wire_spacing,
                length=circuit_wire_length,
                start_x=wire_start_x,
                idx=i + 1,
            )
            wires.append(w)
            self.play(Create(w.line), FadeIn(w.start_circ), FadeIn(w.end_circ), Write(w.label), run_time=0.3)
        
        # PSEUDOCODE (layout stage):
        # - Place H gates on each wire in alternating slots.
        # - Place phase-train blocks between consecutive H gates.
        # - Keep x-bound alignment explicit so later tensor merges are stable.
        # Build circuit layout matching the QFT -> MPO storyline
        gates_and_trains = []
        
        h_x_positions = [element_x_positions[2 * i] for i in range(qubits)]
        pt_x_positions = [element_x_positions[(2 * i) + 1] for i in range(qubits - 1)]

        if len(pt_x_positions) >= 2:
            phase_train_spacing = pt_x_positions[1] - pt_x_positions[0]
            zip_center_x = (pt_x_positions[0] + pt_x_positions[1]) / 2
        else:
            phase_train_spacing = GATE_CENTER_SPACING * 2
            zip_center_x = CIRCUIT_CENTER_X

        PHASE_GATE_WIDTH = GATE_WIDTH
        ZIP_LEFT_INSET = 0.4 * PHASE_GATE_WIDTH
        ZIP_RIGHT_INSET = 0.2 * PHASE_GATE_WIDTH

        zip_x_left = zip_center_x - (phase_train_spacing / 2) + ZIP_LEFT_INSET
        zip_x_right = zip_center_x + (phase_train_spacing / 2) - ZIP_RIGHT_INSET

        phase_palette = [COLOR_BROWN_1, COLOR_BROWN_2, COLOR_BROWN_3]

        def build_phase_targets(control_idx):
            targets = []
            for target_idx in range(control_idx + 1, qubits):
                distance = target_idx - control_idx
                phase_denominator = str(2 ** distance)
                color_idx = min(distance - 1, len(phase_palette) - 1)
                targets.append((wires[target_idx].y_pos, phase_denominator, phase_palette[color_idx]))
            return targets

        h_gates = []
        phase_trains = []

        for control_idx in range(qubits):
            h_gate = QuantumGate("H", h_x_positions[control_idx], wires[control_idx].y_pos, COLOR_H)
            h_gates.append(h_gate)
            gates_and_trains.append(h_gate)

            if control_idx < qubits - 1:
                phase_train = PhaseTrain(
                    pt_x_positions[control_idx],
                    wires[control_idx].y_pos,
                    build_phase_targets(control_idx),
                )
                phase_trains.append(phase_train)
                gates_and_trains.append(phase_train)

        # Compatibility aliases for downstream zip-up / push-down logic
        h1, h2, h3, h4 = h_gates
        pt1, pt2, pt3 = phase_trains
        
        # PSEUDOCODE (initial reveal):
        # - Animate each H as box-grow + label-write.
        # - Animate each phase train as control dot + vertical leg + targets.
        # This establishes the uncompressed baseline circuit.
        # Animate gate creation sequentially
        for item in gates_and_trains:
            if isinstance(item, QuantumGate):
                item.box.save_state()
                item.box.stretch_to_fit_height(0.01)
                self.play(Restore(item.box), Write(item.label), run_time=0.15)
            else: # PhaseTrain
                item.remove(item.cline) # Remove cline so we can create it last
                self.play(GrowFromCenter(item.dot), Create(item.cline), run_time=0.2)
                for g in item.gates:
                    # g[0] is box, g[1] is label
                    g[0].save_state()
                    g[0].stretch_to_fit_height(0.01)
                    self.play(Restore(g[0]), Write(g[1]), run_time=0.15)
                item.add(item.cline)

        # ========================== PHASE 1: ZIP-UP SWEEP ==========================
        # Compress from lower wires to upper wires:
        # local gates -> T_i, then T_i -> (U_iS_i, V_i).
        # ===========================================================================
        self.wait(0.5)
        
        # Factory for wide horizontal tensors (T, US)
        def make_wide_tensor(x_left, x_right, y_pos, label_str, color):
            width = x_right - x_left
            center_x = (x_left + x_right) / 2
            box = RoundedRectangle(corner_radius=0.1, width=width, height=GATE_WIDTH, 
                                   color=COLOR_BLACK, stroke_width=WIRE_THICKNESS, 
                                   fill_color=color, fill_opacity=1).move_to(RIGHT * center_x + UP * y_pos)
            label = MathTex(label_str, color=BLACK, stroke_color=BLACK, stroke_width=1).move_to(box.get_center())
            group = VGroup(box, label)
            group.set_z_index(5)
            return group

        # Factory for vertical tall tensors (used during push-down merge)
        def make_tall_tensor(x_pos, y_top, y_bottom, label_str, color, level_padding=0.4):
            height = abs(y_top - y_bottom) + (2 * level_padding)
            center_y = (y_top + y_bottom) / 2
            box = RoundedRectangle(
                corner_radius=0.1,
                width=GATE_WIDTH,
                height=height,
                color=COLOR_BLACK,
                stroke_width=WIRE_THICKNESS * 1.0,
                fill_color=color,
                fill_opacity=1,
            ).move_to(RIGHT * x_pos + UP * center_y)
            label = MathTex(label_str, color=BLACK, stroke_color=BLACK, stroke_width=1).move_to(box.get_center())
            group = VGroup(box, label)
            group.set_z_index(5)
            return group

        def pulse_highlight_box(target_mob, pulse_count=2, pulse_time=0.3, low_alpha=0.25, high_alpha=0.9):
            highlight = RoundedRectangle(
                corner_radius=0.18,
                width=target_mob.width + 0.4,
                height=target_mob.height + 0.4,
                color=RED,
                stroke_width=5,
                fill_color=RED,
                fill_opacity=0.075,
            ).move_to(target_mob.get_center()).set_z_index(20)

            highlight.set_stroke(opacity=low_alpha)
            self.add(highlight)
            for _ in range(pulse_count):
                self.play(highlight.animate.set_stroke(opacity=high_alpha), run_time=pulse_time)
                self.play(highlight.animate.set_stroke(opacity=low_alpha), run_time=pulse_time)
            self.remove(highlight)

        # Compression-state glossary:
        # current_us: active U_iS_i tensor carried between zip-up steps.
        # current_bond: heavy bond between latest V_i and current_us.
        # wire_objects[i]: mobjects that should translate with wire i during spacing changes.
        # v_nodes / v_map: V tensors captured for later push-down recursion.
        # old_bonds: historical heavy bonds that are shifted/rewired/faded as structure evolves.
        current_us = None
        current_bond = None
        EXPAND_AMOUNT = 1.2  # Space added between wires during SVD split
        HALF = EXPAND_AMOUNT / 2
        
        wire_objects = {i: [] for i in range(qubits)}
        wire_objects[0].extend([h1, pt1.dot])
        wire_objects[1].extend([h2, pt2.dot, pt1.gates[0]])
        wire_objects[2].extend([h3, pt3.dot, pt1.gates[1], pt2.gates[0]])
        wire_objects[3].extend([h4, pt1.gates[2], pt2.gates[1], pt3.gates[0]])
        
        v_nodes = []
        v_map = {}
        old_bonds = []

        # Mapping defining which gates are merged into the Ti tensor for each wire index
        # 3 (n4): merges pt1's and pt2's targets on wire 4
        # 2 (n3): merges the US from the split below + pt1 & pt2 targets on n3
        # 1 (n2): merges the US from the split below + pt1 target on n2 + h2 + pt2's control dot
        zip_target_map = {
            3: lambda us: [pt1.gates[2], pt2.gates[1]],
            2: lambda us: [us, pt1.gates[1], pt2.gates[0]],
            1: lambda us: [us, pt1.gates[0], h2, pt2.dot],
        }

        # List of gates to be explicitly removed from the scene after a merge
        zip_cleanup_map = {
            3: [pt1.gates[2], pt2.gates[1]],
            2: [pt1.gates[1], pt2.gates[0]],
            1: [pt1.gates[0], h2, pt2.dot],
        }

        def remove_zip_sources(w_idx):
            sources = zip_cleanup_map[w_idx]
            if w_idx == 3:
                wire_objects[3].remove(pt1.gates[2])
                wire_objects[3].remove(pt2.gates[1])
                pt1.remove(pt1.gates[2])
                pt2.remove(pt2.gates[1])
            elif w_idx == 2:
                wire_objects[2].remove(pt1.gates[1])
                wire_objects[2].remove(pt2.gates[0])
                pt1.remove(pt1.gates[1])
                pt2.remove(pt2.gates[0])
            elif w_idx == 1:
                wire_objects[1].remove(pt1.gates[0])
                wire_objects[1].remove(h2)
                wire_objects[1].remove(pt2.dot)
                pt1.remove(pt1.gates[0])
                pt2.remove(pt2.dot)
            self.remove(*sources)

        def x_bounds(mob):
            if isinstance(mob, Dot):
                x = mob.get_center()[0]
                return x - (GATE_WIDTH / 2), x + (GATE_WIDTH / 2)
            rect = mob.box if hasattr(mob, "box") else (mob[0] if isinstance(mob, VGroup) and len(mob) > 0 else mob)
            return rect.get_left()[0], rect.get_right()[0]

        def span_x(targets, ignore=None, fallback=(None, None)):
            bounds = [x_bounds(g) for g in targets if g is not None and g is not ignore]
            return (min(l for l, _ in bounds), max(r for _, r in bounds)) if bounds else fallback
        
        for w_idx in [3, 2, 1]:
            # PSEUDOCODE (one zip-up rung):
            # A) Highlight current merge targets on wire w_idx.
            # B) If current_us exists, pre-merge same-wire objects into wire_merge.
            # C) Move wires/bonds to the contracted layout.
            # D) Merge [current_us + wire_merge] (or initial gates) into T_{w_idx+1}.
            # E) Either continue via split (T -> U_iS_i + V_i) or finish pivot logic at w_idx=1.
            target_gates = VGroup(*[g for g in zip_target_map[w_idx](current_us) if g is not None])
            x_A, x_B = span_x(target_gates, ignore=current_us, fallback=(zip_x_left, zip_x_right))
            center_x = (x_A + x_B) / 2

            pulse_highlight_box(target_gates)
            
            wire_merge = None
            if current_us is not None:
                wire_targets = [g for g in target_gates if g is not current_us]
                color_source = next((g[0].get_fill_color() for g in wire_targets if not isinstance(g, Dot)), COLOR_T_US)
                wire_merge = make_wide_tensor(x_A, x_B, wires[w_idx].y_pos, "", color_source)
                wire_merge[1].set_opacity(0)
                # Set z-index for wire_merge to be behind US
                wire_merge.set_z_index(5)

                if w_idx == 1:
                    p_gate = wire_targets[0]   # Ppi/2 gate
                    h_gate = wire_targets[1]   # H gate on wire n2
                    ctrl_dot = wire_targets[2] # Control dot for pt2
                    
                    # Ensure H is visually behind P before moving
                    h_gate.box.set_z_index(p_gate[0].get_z_index() - 1)
                    if hasattr(h_gate, "label"):
                        h_gate.label.set_z_index(p_gate[0].get_z_index() - 1)
                    
                    # 1) Move H behind P and unwrite P's label simultaneously
                    self.play(
                        h_gate.animate.move_to(p_gate[0].get_center()),
                        Unwrite(p_gate[1]),
                        run_time=0.5
                    )
                    
                    # 2) Remove H once covered/merged
                    self.remove(h_gate)
                    
                    # 3) Prepare targets for unified tensor merge
                    copied_p = wire_merge.copy()
                    copied_ctrl = wire_merge.copy()
                    
                    # Bring them to a higher layer for the merge animation
                    # But keep them behind current_us (which will be at z-index 10)
                    ctrl_dot.set_z_index(4)
                    p_gate[0].set_z_index(4)
                    wire_merge.set_z_index(4)
                    if current_us is not None:
                        current_us.set_z_index(10)

                    self.wait(0.2)

                    # 4) Merge the remaining P box and Control dot into the wide T2 tensor
                    self.play(
                        ReplacementTransform(p_gate[0], copied_p[0]),
                        ReplacementTransform(ctrl_dot, copied_ctrl[0]),
                        run_time=0.6
                    )
                    
                    self.remove(copied_p, copied_p[0], copied_p[1])
                    self.remove(copied_ctrl, copied_ctrl[0], copied_ctrl[1])
                    self.remove(p_gate, ctrl_dot)
                    self.add(wire_merge)
                else:
                    copied_wire = [wire_merge.copy() for _ in wire_targets]
                    wire_anim = []
                    for i, g in enumerate(wire_targets):
                        if isinstance(g, Dot):
                            copied_wire[i][0].set_z_index(1)
                            wire_anim.append(ReplacementTransform(g, copied_wire[i][0]))
                        else:
                            g[0].set_z_index(1)
                            g[1].set_z_index(1)
                            # Transform box to the unified wire-merge tensor and "write out" the old label
                            wire_anim.append(ReplacementTransform(g[0], copied_wire[i][0]))
                            wire_anim.append(Unwrite(g[1]))
    
                    self.play(*wire_anim, run_time=0.5)
                    for cw in copied_wire:
                        self.remove(cw, cw[0], cw[1])
                    self.remove(*wire_targets)
                    self.add(wire_merge)
                
                target_gates = VGroup(current_us, wire_merge)

            shrink_anims = []
            new_pt3_line_up = None
            new_pt1_line_up = None
            new_pt2_line_up = None
            if current_us is not None:
                # Shrink -> Wires > w_idx move UP by HALF. Wires <= w_idx move DOWN by HALF.
                # Wires above the current merge point move UP; wires at or below move DOWN
                for wi in range(w_idx + 1, qubits):
                    shrink_anims.append(wires[wi].animate.shift(UP * HALF))
                    for obj in wire_objects[wi]:
                        shrink_anims.append(obj.animate.shift(UP * HALF))
                for wi in range(0, w_idx + 1):
                    shrink_anims.append(wires[wi].animate.shift(DOWN * HALF))
                    for obj in wire_objects[wi]:
                        shrink_anims.append(obj.animate.shift(DOWN * HALF))
                
                # Make sure the wire_merge also move correctly with its wire
                if wire_merge is not None:
                    shrink_anims.append(wire_merge.animate.shift(DOWN * HALF))
                
                for ob in old_bonds:
                    shrink_anims.append(ob.animate.shift(UP * HALF))
                
                if w_idx == 2:
                    new_pt3_line_up = Line(pt3.dot.get_center() + DOWN * HALF, pt3.gates[0][0].get_top() + UP * HALF, stroke_width=WIRE_THICKNESS, color=WHITE).set_z_index(0)
                    shrink_anims.append(ReplacementTransform(pt3.cline, new_pt3_line_up))
                    
                    new_merged_y = wires[w_idx].y_pos - HALF
                    new_pt1_line_up = Line(pt1.dot.get_center() + DOWN * HALF, RIGHT * pt1.dot.get_center()[0] + UP * (new_merged_y + 0.4), stroke_width=WIRE_THICKNESS, color=WHITE).set_z_index(0)
                    new_pt2_line_up = Line(pt2.dot.get_center() + DOWN * HALF, RIGHT * pt2.dot.get_center()[0] + UP * (new_merged_y + 0.4), stroke_width=WIRE_THICKNESS, color=WHITE).set_z_index(0)
                    shrink_anims.extend([ReplacementTransform(pt1.cline, new_pt1_line_up), ReplacementTransform(pt2.cline, new_pt2_line_up)])
                elif w_idx == 1:
                    shrink_anims.append(pt3.cline.animate.shift(UP * HALF))
                    # pt2.cline merged into T2 completely, so we remove the dangling leg!
                    shrink_anims.append(FadeOut(pt2.cline))
                    
                    new_merged_y = wires[w_idx].y_pos - HALF
                    new_pt1_line_up = Line(pt1.dot.get_center() + DOWN * HALF, RIGHT * pt1.dot.get_center()[0] + UP * (new_merged_y + 0.4), stroke_width=WIRE_THICKNESS, color=WHITE).set_z_index(0)
                    shrink_anims.append(ReplacementTransform(pt1.cline, new_pt1_line_up))

            merged_y = wires[w_idx].y_pos if current_us is None else wires[w_idx].y_pos - HALF
            merged_tensor = make_wide_tensor(x_A, x_B, merged_y, f"T_{w_idx+1}", COLOR_T_US)
            
            source_targets = list(target_gates)
            us_moving_label = None
            us_moving_box = None
            if current_us is not None:
                # Pause at the actual U_iS_i + wire-merge -> T_i step.
                self.wait(0.4)
                merged_tensor[1].set_opacity(0)
                us_moving_label = current_us[1].copy().set_z_index(11)
                us_moving_box = current_us[0]
                us_moving_label.add_updater(lambda m: m.move_to(us_moving_box.get_center()))
                current_us[1].set_opacity(0)
                # Ensure US is on top during the merge
                current_us.set_z_index(10)
                # wire_merge.set_z_index(5)
                self.add(us_moving_label)
            else:
                merged_tensor[1].set_opacity(0)

            copied_targets = [merged_tensor.copy() for _ in target_gates]
            anim_group = list(shrink_anims)
            if current_us is not None:
                label_source_idx = None
            else:
                label_source_idx = next((i for i, g in enumerate(target_gates) if not isinstance(g, Dot)), None)
            
            for i, g in enumerate(target_gates):
                if isinstance(g, Dot):
                    # For Dot, transform the dot itself into the target box
                    anim_group.append(ReplacementTransform(g, copied_targets[i][0]))
                else:
                    # For VGroup (TensorNode or wire_merge), transform box to box
                    # This prevents ValueError: zip() argument 2 is shorter than argument 1
                    # which happens when MathTex labels have different submobject counts.
                    if label_source_idx is None or i != label_source_idx:
                        copied_targets[i][1].set_opacity(0)
                    
                    if current_us is not None and i == 0:
                        # Transforming current_us (i=0)
                        anim_group.append(ReplacementTransform(g[0], copied_targets[i][0]))
                        anim_group.append(Unwrite(us_moving_label))
                    elif current_us is None:
                        # Transforming initial gates (w_idx=3)
                        anim_group.append(ReplacementTransform(g[0], copied_targets[i][0]))
                        anim_group.append(Unwrite(g[1]))
                    else:
                        # Transforming wire_merge (i=1)
                        # No label unwrite needed as it was empty/hidden
                        anim_group.append(ReplacementTransform(g[0], copied_targets[i][0]))
            
            if current_bond is not None:
                new_old_bond = Line(v_nodes[-1].box.get_top() + UP * HALF, merged_tensor[0].get_bottom(), color=COLOR_LIGHT_BLUE, stroke_width=BOND_LIGHT).set_z_index(-1)
                anim_group.append(ReplacementTransform(current_bond, new_old_bond))
                current_bond = new_old_bond
                old_bonds.append(current_bond)
                
            self.play(*anim_group, run_time=0.5)
            self.remove(*source_targets)
            if us_moving_label is not None:
                us_moving_label.clear_updaters()
                self.remove(us_moving_label)
            for ct in copied_targets:
                self.remove(ct, ct[0], ct[1])
            self.add(merged_tensor)
            merged_tensor[1].set_opacity(1)
            merged_tensor[1].set_z_index(6)
            self.play(Write(merged_tensor[1]), run_time=0.22)
            
            if current_us is not None:
                if w_idx == 2:
                    pt3.cline = new_pt3_line_up
                    pt1.cline = new_pt1_line_up
                    pt2.cline = new_pt2_line_up
                elif w_idx == 1:
                    pt1.cline = new_pt1_line_up
                for wi in range(w_idx + 1, qubits):
                    wires[wi].y_pos += HALF
                for wi in range(0, w_idx + 1):
                    wires[wi].y_pos -= HALF
            
            remove_zip_sources(w_idx)

            if w_idx == 1:
                # PSEUDOCODE (top pivot in zip-up):
                # - Build U_1 from H_1 + control dot.
                # - Convert merged T to centered T_2 with updated thin bond.
                # - Recenter lower structures to maintain overall frame balance.
                # U1 Creation
                u1_targets = VGroup(h1, pt1.dot)
                pulse_highlight_box(u1_targets)
                control_x = pt1.dot.get_center()[0]
                
                u1_initial = TensorNode("U_1", control_x, wires[0].y_pos, "U").set_z_index(5)
                copied_u1s = [u1_initial.copy() for _ in u1_targets]
                
                self.play(
                    ReplacementTransform(h1, copied_u1s[0]),
                    ReplacementTransform(pt1.dot, copied_u1s[1][0]),
                    run_time=0.6
                )
                for cu in copied_u1s:
                    self.remove(cu, cu[0], cu[1])
                self.add(u1_initial)
                
                # Reshape and Center
                t_center_sq = TensorNode("T_2", control_x, merged_y, "T").set_z_index(5)
                t_center_sq.box.set_fill(COLOR_T_US)
                
                u1_final = TensorNode("U_1", control_x, wires[0].y_pos, "U").set_z_index(5)
                # Keep bond thin and white to maintain original dimension!
                final_bond = Line(u1_final.box.get_bottom(), t_center_sq.box.get_top(), stroke_width=WIRE_THICKNESS, color=WHITE).set_z_index(0)

                anchor_x = u1_final.box.get_center()[0]
                ref_x = current_us[0].get_center()[0] if current_us is not None else anchor_x
                dx_lower = anchor_x - ref_x
                object_dx_map = {}
                lower_objs = []
                for wi in range(2, qubits):
                    lower_objs.extend(wire_objects[wi])
                lower_objs.extend([pt3.cline, current_us, current_bond, *old_bonds])

                extra_right_dx = 0.0
                if h3 in self.mobjects:
                    tt_right_after = t_center_sq.box.get_right()[0]
                    h3_left_after = h3.box.get_left()[0] + dx_lower
                    extra_right_dx = (tt_right_after + MIN_GATE_GAP) - h3_left_after

                right_block = [h3, pt3.dot, pt3.gates[0], h4, pt3.cline]
                scene_family_ids = {id(m) for root in self.mobjects for m in root.get_family()}
                seen_ids = set()
                for obj in lower_objs:
                    if obj is None or id(obj) not in scene_family_ids or id(obj) in seen_ids:
                        continue
                    total_dx = dx_lower + (extra_right_dx if obj in right_block else 0.0)
                    object_dx_map[id(obj)] = total_dx
                    seen_ids.add(id(obj))

                wire_family_ids = {id(m) for w in wires for m in w.get_family()}
                movable_objs = [obj for obj in self.mobjects if id(obj) not in wire_family_ids]
                replace_sources = {id(merged_tensor), id(u1_initial), id(pt1.cline)}

                predicted_bounds = []
                for obj in movable_objs:
                    if id(obj) in replace_sources:
                        continue
                    dx_obj = object_dx_map.get(id(obj), 0.0)
                    predicted_bounds.append((obj.get_left()[0] + dx_obj, obj.get_right()[0] + dx_obj))
                predicted_bounds.extend([
                    (t_center_sq.box.get_left()[0], t_center_sq.box.get_right()[0]),
                    (u1_final.box.get_left()[0], u1_final.box.get_right()[0]),
                    (final_bond.get_left()[0], final_bond.get_right()[0]),
                ])

                recenter_dx = 0.0
                if predicted_bounds:
                    predicted_left = min(left for left, _ in predicted_bounds)
                    predicted_right = max(right for _, right in predicted_bounds)
                    recenter_dx = CIRCUIT_CENTER_X - ((predicted_left + predicted_right) / 2)

                t_center_sq.shift(RIGHT * recenter_dx)
                u1_final.shift(RIGHT * recenter_dx)
                final_bond.shift(RIGHT * recenter_dx)

                shift_anims = []
                for obj in movable_objs:
                    if id(obj) in replace_sources:
                        continue
                    total_dx = object_dx_map.get(id(obj), 0.0) + recenter_dx
                    if abs(total_dx) > 1e-6:
                        shift_anims.append(obj.animate.shift(RIGHT * total_dx))
                
                self.play(
                    ReplacementTransform(merged_tensor, t_center_sq),
                    ReplacementTransform(u1_initial, u1_final),
                    ReplacementTransform(pt1.cline, final_bond),
                    *shift_anims,
                    run_time=0.8
                )
                self.remove(t_center_sq.box, t_center_sq.label, u1_final.box, u1_final.label)
                self.remove(u1_initial, u1_initial[0], u1_initial[1])
                self.remove(merged_tensor, merged_tensor[0], merged_tensor[1])
                self.add(t_center_sq, u1_final)
                
                t_center = t_center_sq
                break

            # ZIP-UP continuation:
            # Expand wire spacing and split merged tensor into U_iS_i + V_i.
            shift_down_anims = []
            for wi in range(w_idx, qubits):
                shift_down_anims.append(wires[wi].animate.shift(DOWN * HALF))
                for obj in wire_objects[wi]:
                    shift_down_anims.append(obj.animate.shift(DOWN * HALF))
            for wi in range(0, w_idx):
                shift_down_anims.append(wires[wi].animate.shift(UP * HALF))
                for obj in wire_objects[wi]:
                    shift_down_anims.append(obj.animate.shift(UP * HALF))
            
            for ob in old_bonds[:-1] if current_bond is not None else old_bonds:
                shift_down_anims.append(ob.animate.shift(DOWN * HALF))
            
            new_pt3_line_down = None
            if w_idx == 3:
                new_pt3_line_down = Line(pt3.dot.get_center() + UP * HALF, pt3.gates[0][0].get_top() + DOWN * HALF, stroke_width=WIRE_THICKNESS, color=WHITE).set_z_index(0)
                shift_down_anims.append(ReplacementTransform(pt3.cline, new_pt3_line_down))
            else:
                shift_down_anims.append(pt3.cline.animate.shift(DOWN * HALF))

            v_y = wires[w_idx].y_pos - HALF
            v_node = TensorNode(f"V_{w_idx+1}", center_x, v_y, "V").set_z_index(5)
            
            # Place the US node midway between the original wires being separated
            us_y = (v_y + wires[w_idx-1].y_pos + HALF) / 2
            us_node = make_wide_tensor(x_A, x_B, us_y, f"U_{w_idx+1}S_{w_idx+1}", COLOR_T_US)
            
            fat_bond = Line(v_node.box.get_top(), us_node[0].get_bottom(), color=COLOR_ORANGE_RED, stroke_width=BOND_HEAVY).set_z_index(-1)
            
            pt1_x = pt1.dot.get_center()[0]
            pt2_x = pt2.dot.get_center()[0]
            new_pt1_cline = Line(
                pt1.dot.get_center() + UP * HALF,
                RIGHT * pt1_x + UP * (us_y + 0.4),
                stroke_width=WIRE_THICKNESS,
                color=WHITE,
            ).set_z_index(0)
            new_pt2_cline = Line(
                pt2.dot.get_center() + UP * HALF,
                RIGHT * pt2_x + UP * (us_y + 0.4),
                stroke_width=WIRE_THICKNESS,
                color=WHITE,
            ).set_z_index(0)

            bond_anims = []
            if current_bond is not None:
                new_old_bond2 = Line(v_nodes[-1].box.get_top() + DOWN * HALF, v_node.box.get_bottom(), color=COLOR_LIGHT_BLUE, stroke_width=BOND_LIGHT).set_z_index(-1)
                bond_anims.append(ReplacementTransform(current_bond, new_old_bond2))
                old_bonds[-1] = new_old_bond2

            # Hold on T_i briefly before showing the U_i S_i / V_i split.
            self.wait(0.4)

            self.play(
                ReplacementTransform(merged_tensor, us_node),
                Create(v_node.box),
                Write(v_node.label),
                Create(fat_bond),
                ReplacementTransform(pt1.cline, new_pt1_cline),
                ReplacementTransform(pt2.cline, new_pt2_cline),
                *bond_anims,
                *shift_down_anims,
                run_time=0.8
            )
            self.remove(v_node.box, v_node.label, us_node[0], us_node[1])
            self.remove(merged_tensor, merged_tensor[0], merged_tensor[1])
            self.add(v_node, us_node)
            
            pt1.cline = new_pt1_cline
            pt2.cline = new_pt2_cline
            if w_idx == 3: pt3.cline = new_pt3_line_down
            
            current_us = us_node
            current_bond = fat_bond
            
            wire_objects[w_idx].append(v_node)
            v_nodes.append(v_node)
            v_map[w_idx + 1] = v_node
            
            for wi in range(w_idx, qubits): wires[wi].y_pos -= HALF
            for wi in range(0, w_idx): wires[wi].y_pos += HALF

            self.wait(0.4)
            thin_bond = Line(v_node.box.get_top(), us_node[0].get_bottom(), color=COLOR_LIGHT_BLUE, stroke_width=BOND_LIGHT).set_z_index(-1)
            self.play(ReplacementTransform(fat_bond, thin_bond), run_time=0.8)
            current_bond = thin_bond

        self.wait(1)
        
        # ======================== PHASE 2: PUSH-DOWN SWEEP ========================
        # Recursive orthogonalization from top wires toward bottom wires.
        # ==========================================================================
        down_x = u1_final.box.get_center()[0] if 'u1_final' in locals() else t_center.box.get_center()[0]

        # Recursive function to handle the downward decomposition
        def compress_down_recursive(train_idx, t_current, previous_u=None, blue_bond_to_current_t=None):
            # PSEUDOCODE (recursive unit):
            # 1) Merge (t_current, V_{train_idx+1}) into a tall temporary tensor.
            # 2) Split into (U_train_idx, T_{train_idx+1}).
            # 3) Reconnect compressed blue bonds from previous U (if any).
            # 4) Recurse until no next V exists.
            axis_x = t_current.box.get_center()[0]

            if train_idx >= qubits:
                return t_current, previous_u, blue_bond_to_current_t

            v_next = v_map.get(train_idx + 1)
            if v_next is None:
                return t_current, previous_u, blue_bond_to_current_t

            pulse_highlight_box(VGroup(t_current, v_next))

            merged_pair = make_tall_tensor(
                axis_x,
                wires[train_idx - 1].y_pos,
                wires[train_idx].y_pos,
                "",
                COLOR_T_US,
                level_padding=0.55,
            )
            merged_pair[1].set_opacity(0)

            merge_anims = [
                ReplacementTransform(t_current.box, merged_pair[0]),
                Unwrite(t_current.label),
                Unwrite(v_next.label),
                FadeOut(v_next.box),
            ]

            if train_idx == 2 and current_bond is not None:
                merge_anims.append(FadeOut(current_bond))
            if train_idx > 2 and old_bonds:
                merge_anims.append(FadeOut(old_bonds[0]))

            blue_bond_to_merged = None
            if previous_u is not None and blue_bond_to_current_t is not None:
                blue_bond_to_merged = Line(previous_u.box.get_bottom(), merged_pair[0].get_top(), color=COLOR_LIGHT_BLUE, stroke_width=BOND_LIGHT).set_z_index(0)
                merge_anims.append(ReplacementTransform(blue_bond_to_current_t, blue_bond_to_merged))

            self.play(*merge_anims, run_time=0.42)
            self.remove(merged_pair[0], merged_pair[1])
            self.add(merged_pair)
            self.wait(0.2)

            u_node = TensorNode(f"U_{train_idx}", axis_x, wires[train_idx - 1].y_pos, "U")
            t_next = TensorNode(f"T_{train_idx+1}", axis_x, wires[train_idx].y_pos, "T")
            t_next.box.set_fill(COLOR_T_US)
            vertical_blue_bond = Line(u_node.box.get_bottom(), t_next.box.get_top(), color=COLOR_LIGHT_BLUE, stroke_width=BOND_LIGHT).set_z_index(0)

            split_anims = [
                ReplacementTransform(merged_pair[0], u_node.box),
                FadeOut(merged_pair[1]),
                Write(u_node.label),
                Create(t_next.box),
                Write(t_next.label),
                Create(vertical_blue_bond),
            ]

            if previous_u is not None and blue_bond_to_merged is not None:
                split_anims.append(
                    ReplacementTransform(
                        blue_bond_to_merged,
                        Line(previous_u.box.get_bottom(), u_node.box.get_top(), color=COLOR_LIGHT_BLUE, stroke_width=BOND_LIGHT).set_z_index(0),
                    )
                )

            self.play(*split_anims, run_time=0.8)
            self.remove(u_node.box, u_node.label, t_next.box, t_next.label)
            self.remove(merged_pair, merged_pair[0], merged_pair[1])
            self.add(u_node, t_next, vertical_blue_bond)

            return compress_down_recursive(train_idx + 1, t_next, u_node, vertical_blue_bond)

        # ================= PHASE 3: CUSTOM RIGHT-EDGE CHOREOGRAPHY =================
        # Explicitly handle bottom/right merges (n4 then n3), including control-leg
        # rerouting and staged recenter shifts that are more complex than the generic
        # recursive path.
        # ===========================================================================
        # Handles the rightward horizontal merge on the bottom wire (n4)
        bottom_t, upper_u, upper_to_bottom_bond = compress_down_recursive(2, t_center)
        phase_on_bottom = pt3.gates[0] if len(pt3.gates) > 0 and pt3.gates[0] in self.mobjects else None
        upper_h = h3 if h3 in self.mobjects else None
        upper_control = pt3.dot if pt3.dot in self.mobjects else None
        right_h = h4 if h4 in self.mobjects else None

        if phase_on_bottom is None:
            return

        pulse_highlight_box(VGroup(bottom_t, phase_on_bottom))

        t4_left = min(bottom_t.box.get_left()[0], phase_on_bottom[0].get_left()[0])
        t4_right = max(bottom_t.box.get_right()[0], phase_on_bottom[0].get_right()[0])
        t4_phase_merge = make_wide_tensor(t4_left, t4_right, wires[3].y_pos, "T_4", COLOR_T_US)
        t4_phase_copy = t4_phase_merge.copy()

        self.play(
            ReplacementTransform(bottom_t, t4_phase_merge),
            ReplacementTransform(phase_on_bottom, t4_phase_copy),
            run_time=0.55,
        )
        self.remove(t4_phase_copy, t4_phase_copy[0], t4_phase_copy[1])
        self.remove(bottom_t, phase_on_bottom)
        if phase_on_bottom in wire_objects[3]:
            wire_objects[3].remove(phase_on_bottom)
        if phase_on_bottom in pt3.gates:
            pt3.remove(phase_on_bottom)

        control_to_t_leg = None
        if upper_control is not None and pt3.cline in self.mobjects:
            cx = upper_control.get_center()[0]
            control_to_t_leg = Line(
                upper_control.get_center(),
                RIGHT * cx + UP * t4_phase_merge[0].get_top()[1],
                stroke_width=WIRE_THICKNESS,
                color=WHITE,
            ).set_z_index(0)
            self.play(ReplacementTransform(pt3.cline, control_to_t_leg), run_time=0.2)
            pt3.cline = control_to_t_leg

        expand = 1.2
        us4_y = (wires[2].y_pos + wires[3].y_pos) / 2
        us4 = make_wide_tensor(
            t4_phase_merge[0].get_left()[0],
            t4_phase_merge[0].get_right()[0],
            us4_y,
            "U_4S_4",
            COLOR_T_US,
        )

        v4_x = us4[0].get_center()[0]
        v4_y = wires[3].y_pos - (expand / 2)
        v4_new = TensorNode("V_4", v4_x, v4_y, "V")
        red_34 = Line(us4[0].get_bottom(), v4_new.box.get_top(), color=COLOR_ORANGE_RED, stroke_width=BOND_HEAVY).set_z_index(0)

        u3_x = upper_u.box.get_center()[0]
        blue_u3_us4 = Line(
            RIGHT * u3_x + UP * (upper_u.box.get_bottom()[1] + (expand / 2)),
            RIGHT * u3_x + UP * us4[0].get_top()[1],
            color=COLOR_LIGHT_BLUE,
            stroke_width=BOND_LIGHT,
        ).set_z_index(0)

        control_leg_to_us4 = None
        if control_to_t_leg is not None:
            cx = upper_control.get_center()[0]
            control_leg_to_us4 = Line(
                RIGHT * cx + UP * (upper_control.get_center()[1] + (expand / 2)),
                RIGHT * cx + UP * us4[0].get_top()[1],
                stroke_width=WIRE_THICKNESS,
                color=WHITE,
            ).set_z_index(0)

        expand_down_objs = [obj for obj in [right_h] if obj is not None and obj in self.mobjects]
        recenter_objs = [
            mob for mob in list(self.mobjects)
            if mob not in expand_down_objs and mob is not wires[3] and mob is not upper_to_bottom_bond and mob is not control_to_t_leg and mob is not t4_phase_merge
        ]

        self.play(
            wires[3].animate.shift(DOWN * (expand / 2)),
            *[obj.animate.shift(DOWN * (expand / 2)) for obj in expand_down_objs],
            *[mob.animate.shift(UP * (expand / 2)) for mob in recenter_objs],
            ReplacementTransform(t4_phase_merge, us4),
            Create(v4_new.box),
            Write(v4_new.label),
            Create(red_34),
            ReplacementTransform(upper_to_bottom_bond, blue_u3_us4),
            *([ReplacementTransform(control_to_t_leg, control_leg_to_us4)] if control_to_t_leg is not None else []),
            run_time=0.6,
        )
        self.remove(v4_new.box, v4_new.label)
        self.add(v4_new)
        upper_to_bottom_bond = blue_u3_us4
        if control_leg_to_us4 is not None:
            control_to_t_leg = control_leg_to_us4
        wires[2].y_pos += expand / 2
        wires[3].y_pos -= expand / 2
        
        self.wait(0.2)
        thin_red_34 = Line(us4[0].get_bottom(), v4_new.box.get_top(), color=COLOR_LIGHT_BLUE, stroke_width=BOND_LIGHT).set_z_index(0)
        self.play(ReplacementTransform(red_34, thin_red_34), run_time=0.6)
        red_34 = thin_red_34

        if upper_h is None or upper_control is None:
            return

        # PSEUDOCODE (n3 pre-merge):
        # Merge [U_3, H_3, control-dot] into n3_merge, then perform:
        # [n3_merge + U_4S_4] -> T_3, while collapsing temporary bonds and
        # restoring wire spacing.
        targets_n3 = VGroup(upper_u, upper_h, upper_control, us4)
        self.wait(0.25)
        pulse_highlight_box(targets_n3)

        n3_merge = make_wide_tensor(us4[0].get_left()[0], us4[0].get_right()[0], wires[2].y_pos, "", upper_u.box.get_fill_color())
        n3_merge[1].set_opacity(0)
        # 1) Move H behind U3 and unwrite U3 label
        upper_h.box.set_z_index(upper_u[0].get_z_index() - 1)
        if hasattr(upper_h, "label"):
            upper_h.label.set_z_index(upper_u[0].get_z_index() - 1)
            
        self.play(
            upper_h.animate.move_to(upper_u.get_center()),
            Unwrite(upper_u[1]),
            run_time=0.45
        )
        self.remove(upper_h)

        # 2) Merge U3 and Control into a yellow intermediate n3_merge
        n3_merge = make_wide_tensor(us4[0].get_left()[0], us4[0].get_right()[0], wires[2].y_pos, "", upper_u.box.get_fill_color())
        n3_merge[1].set_opacity(0)
        n3_merge.set_z_index(4)
        
        n3_tmp_copies = [n3_merge.copy(), n3_merge.copy()]
        upper_u[0].set_z_index(4)
        upper_control.set_z_index(4)

        self.wait(0.2)
        self.play(
            ReplacementTransform(upper_u[0], n3_tmp_copies[0][0]),
            ReplacementTransform(upper_control, n3_tmp_copies[1][0]),
            run_time=0.5,
        )
        for cp in n3_tmp_copies: self.remove(cp, cp[0], cp[1])
        self.remove(upper_u, upper_control)
        self.add(n3_merge)
        self.wait(0.4)

        # 3) Merge n3_merge and US4 into final light-blue T3
        t3_wide = make_wide_tensor(us4[0].get_left()[0], us4[0].get_right()[0], wires[2].y_pos - (expand / 2), "T_3", COLOR_T_US)
        t3_wide[1].set_opacity(0)
        t3_from_n3 = t3_wide.copy()
        t3_from_us4 = t3_wide.copy()

        us4_moving_label = us4[1].copy().set_z_index(11)
        us4_moving_box = us4[0]
        us4_moving_label.add_updater(lambda m: m.move_to(us4_moving_box.get_center()))
        us4[1].set_opacity(0)
        us4.set_z_index(10)
        
        self.add(us4_moving_label)
        restore_n4_objs = [obj for obj in [v4_new, right_h] if obj is not None and obj in self.mobjects]
        recenter_exclusions = {
            id(wires[3]), id(red_34), id(us4), id(n3_merge),
            id(upper_to_bottom_bond), id(control_to_t_leg), id(pt3.cline), id(us4_moving_label)
        }
        restore_recenter_objs = [mob for mob in list(self.mobjects) if mob not in restore_n4_objs and id(mob) not in recenter_exclusions]
        
        red_t3_v4_restored = Line(
            RIGHT * t3_wide[0].get_center()[0] + UP * t3_wide[0].get_bottom()[1],
            RIGHT * v4_new.box.get_center()[0] + UP * (v4_new.box.get_top()[1] + (expand / 2)),
            color=COLOR_LIGHT_BLUE, stroke_width=BOND_LIGHT,
        ).set_z_index(0)

        control_bond_collapse_anims = []
        collapsed_control_bond = None
        collapsed_pt3_cline = None
        if control_to_t_leg is not None and control_to_t_leg in self.mobjects:
            cx, cy = control_to_t_leg.get_center()[0], t3_wide[0].get_center()[1]
            collapsed_control_bond = Line(RIGHT * cx + UP * (cy+0.03), RIGHT * cx + UP * (cy-0.03), color=WHITE, stroke_width=WIRE_THICKNESS).set_z_index(0)
            control_bond_collapse_anims = [ReplacementTransform(control_to_t_leg, collapsed_control_bond)]
        elif pt3.cline in self.mobjects:
            cx, cy = pt3.cline.get_center()[0], t3_wide[0].get_center()[1]
            collapsed_pt3_cline = Line(RIGHT * cx + UP * (cy+0.03), RIGHT * cx + UP * (cy-0.03), color=WHITE, stroke_width=WIRE_THICKNESS).set_z_index(0)
            control_bond_collapse_anims = [ReplacementTransform(pt3.cline, collapsed_pt3_cline)]

        bx, by = upper_to_bottom_bond.get_center()[0], t3_wide[0].get_center()[1]
        collapsed_upper_bond = Line(RIGHT * bx + UP * (by+0.03), RIGHT * bx + UP * (by-0.03), color=COLOR_LIGHT_BLUE, stroke_width=BOND_LIGHT).set_z_index(0)

        self.play(
            ReplacementTransform(n3_merge[0], t3_from_n3[0]),
            ReplacementTransform(us4[0], t3_from_us4[0]),
            Unwrite(us4_moving_label),
            ReplacementTransform(upper_to_bottom_bond, collapsed_upper_bond),
            wires[3].animate.shift(UP * (expand / 2)),
            *[obj.animate.shift(UP * (expand / 2)) for obj in restore_n4_objs],
            *[mob.animate.shift(DOWN * (expand / 2)) for mob in restore_recenter_objs],
            ReplacementTransform(red_34, red_t3_v4_restored),
            *control_bond_collapse_anims,
            run_time=0.6,
        )
        us4_moving_label.clear_updaters()
        self.remove(us4_moving_label, n3_merge, us4, collapsed_upper_bond)
        self.remove(t3_from_n3, t3_from_n3[0], t3_from_n3[1], t3_from_us4, t3_from_us4[0], t3_from_us4[1])
        if collapsed_control_bond: self.remove(collapsed_control_bond)
        if collapsed_pt3_cline: self.remove(collapsed_pt3_cline)
        if collapsed_control_bond: self.remove(collapsed_control_bond)
        if collapsed_pt3_cline: self.remove(collapsed_pt3_cline)
        
        self.add(t3_wide)
        t3_wide[1].set_opacity(1)
        t3_wide[1].set_z_index(6)
        self.play(Write(t3_wide[1]), run_time=0.22)
        red_t3_v4 = red_t3_v4_restored
        wires[2].y_pos -= expand / 2
        wires[3].y_pos += expand / 2

        # Normalize temporary wide T_3 into canonical TensorNode("T_3")
        t3_new = TensorNode("T_3", down_x, wires[2].y_pos, "T")
        t3_new.box.set_fill(COLOR_T_US)
        v4_shift_dx = down_x - t3_wide[0].get_center()[0]

        h4_move_dx = 0.0
        h4_move_dy = 0.0
        if right_h is not None and right_h in self.mobjects:
            h_target_x = down_x + (GATE_WIDTH / 2) + MIN_GATE_GAP + (GATE_WIDTH / 2)
            h_target_y = wires[3].y_pos
            h4_move_dx = h_target_x - right_h.get_center()[0]
            h4_move_dy = h_target_y - right_h.get_center()[1]

        wire_family_ids = {id(m) for w in wires for m in w.get_family()}
        movable_objs = [mob for mob in self.mobjects if id(mob) not in wire_family_ids]
        replace_sources = {id(t3_wide)}

        predicted_bounds = []
        for obj in movable_objs:
            if id(obj) in replace_sources:
                continue
            dx_obj = 0.0
            if obj is v4_new or obj is red_t3_v4:
                dx_obj = v4_shift_dx
            elif obj is right_h:
                dx_obj = h4_move_dx
            predicted_bounds.append((obj.get_left()[0] + dx_obj, obj.get_right()[0] + dx_obj))
        predicted_bounds.append((t3_new.box.get_left()[0], t3_new.box.get_right()[0]))

        recenter_dx = 0.0
        if predicted_bounds:
            predicted_left = min(left for left, _ in predicted_bounds)
            predicted_right = max(right for _, right in predicted_bounds)
            recenter_dx = CIRCUIT_CENTER_X - ((predicted_left + predicted_right) / 2)

        t3_new.shift(RIGHT * recenter_dx)

        recenter_anims = []
        for obj in movable_objs:
            if id(obj) in replace_sources or obj is v4_new or obj is red_t3_v4 or obj is right_h:
                continue
            if abs(recenter_dx) > 1e-6:
                recenter_anims.append(obj.animate.shift(RIGHT * recenter_dx))

        self.play(
            ReplacementTransform(t3_wide[0], t3_new.box),
            ReplacementTransform(t3_wide[1], t3_new.label),
            v4_new.animate.shift(RIGHT * (v4_shift_dx + recenter_dx)),
            red_t3_v4.animate.shift(RIGHT * (v4_shift_dx + recenter_dx)),
            *([right_h.animate.shift(RIGHT * (h4_move_dx + recenter_dx) + UP * h4_move_dy)] if right_h is not None and right_h in self.mobjects else []),
            *recenter_anims,
            run_time=0.45,
        )
        self.remove(t3_wide, t3_wide[0], t3_wide[1])
        self.remove(t3_new.box, t3_new.label)
        self.add(t3_new)

        # ================= PHASE 4: COMPLETE FINAL MPO =================
        # Skip the recursive T3 -> T4 step. Instead, merge the final H gate
        # (h4) directly into V4 and present as a 4-site MPO.
        # ===============================================================
        if right_h is not None:
            pulse_highlight_box(VGroup(v4_new, right_h))
            
            # Temporary copy for the transformation
            h_into_v4 = v4_new.copy()
            h_into_v4[1].set_opacity(0)  # Hide label during transform
            
            self.play(
                ReplacementTransform(right_h, h_into_v4),
                run_time=0.45,
            )
            # Remove the temporary objects; v4_new remains as the final Site 4 tensor
            self.remove(h_into_v4, h_into_v4[0], h_into_v4[1])
            self.remove(right_h)
            if right_h in wire_objects[3]:
                wire_objects[3].remove(right_h)

        self.wait(0.6)

        # ==================== PHASE 5: FINAL MPO DISPLAY ====================
        # Recenter the final U1, U2, T3, V4 MPO and shrink wires to MPO legs.
        # ====================================================================
        
        # Gather all non-wire objects (tensors, labels, bonds)
        wire_family_ids = {id(m) for w in wires for m in w.get_family()}
        mpo_objs = [mob for mob in self.mobjects if id(mob) not in wire_family_ids]
        
        # Group tensors specifically to calculate centering
        tensors_vg = VGroup(*[mob for mob in mpo_objs if hasattr(mob, "box") or isinstance(mob, TensorNode)])
        recenter_dx = CIRCUIT_CENTER_X - tensors_vg.get_center()[0]
        
        # Wire leg settings
        center_x = CIRCUIT_CENTER_X
        half_wire_final = 1.0  # Even longer MPO leg length as requested
        left_x = center_x - half_wire_final
        right_x = center_x + half_wire_final

        final_anims = []
        # Shift all Mpo objects to center
        if abs(recenter_dx) > 1e-6:
            final_anims.extend([mob.animate.shift(RIGHT * recenter_dx) for mob in mpo_objs])

        # Shorten wires and move labels
        for w in wires:
            y = w.line.get_center()[1]
            final_anims.append(w.line.animate.put_start_and_end_on(np.array([left_x, y, 0]), np.array([right_x, y, 0])))
            final_anims.append(w.start_circ.animate.move_to(np.array([left_x, y, 0])))
            final_anims.append(w.end_circ.animate.move_to(np.array([right_x, y, 0])))
            final_anims.append(w.label.animate.move_to(np.array([left_x - 0.4, y, 0])))

        self.play(*final_anims, run_time=1.0)
        self.wait(2.0)
