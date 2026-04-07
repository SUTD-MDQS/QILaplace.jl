from manim import *

config.background_color = "#222831"

# Shared visual constants
GATE_WIDTH = 0.50
WIRE_THICKNESS = 3.5
MIN_GATE_GAP = 0.12
GATE_CENTER_SPACING = GATE_WIDTH + MIN_GATE_GAP
GATE_LABEL_SCALE = 0.5
WIRE_LABEL_SCALE = 0.56
HATCH_LINES_SQUARE_GATE = 14

WIRE_SPACING_MAIN_COPY = 0.62
WIRE_SPACING_PAIR = 0.88
CIRCUIT_CENTER_X = 0.0

COLOR_BLACK = BLACK
COLOR_H = "#3CB371"
COLOR_HD = COLOR_H
COLOR_BROWN_1 = "#D98E43"
COLOR_BROWN_2 = "#A47243"
COLOR_BROWN_3 = "#6A4329"


def set_text_defaults() -> None:
    Text.set_default(color="#DDDDDD")
    MathTex.set_default(color="#DDDDDD")
    Tex.set_default(color="#DDDDDD")


class QuantumWire(VGroup):
    def __init__(
        self,
        y_pos,
        length,
        start_x,
        idx,
        is_copy=False,
        display_idx=None,
        label_scale=WIRE_LABEL_SCALE,
    ):
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
        if display_idx is None:
            display_idx = (idx // 2) + 1
        self.label = MathTex(f"|n_{{{display_idx}{prime}}}\\rangle").scale(label_scale).next_to(self.start_circ, LEFT, buff=0.10)

        self.add(self.line, self.start_circ, self.end_circ, self.label)
        self.set_z_index(-2)


class QuantumGate(VGroup):
    def __init__(self, label_str, x_pos, y_pos, color, label_scale=GATE_LABEL_SCALE):
        super().__init__()
        side = GATE_WIDTH
        self.box = RoundedRectangle(
            corner_radius=0.06,
            width=side,
            height=side,
            color=COLOR_BLACK,
            stroke_width=WIRE_THICKNESS * side,
            fill_color=color,
            fill_opacity=1,
        ).move_to(RIGHT * x_pos + UP * y_pos)
        self.label = MathTex(label_str, color=BLACK, stroke_color=BLACK, stroke_width=0.8).scale(label_scale).move_to(self.box.get_center())
        self.add(self.box, self.label)


class HatchedGate(VGroup):
    def __init__(self, label_str, x_pos, y_pos, color, label_scale=GATE_LABEL_SCALE):
        super().__init__()
        side = GATE_WIDTH
        self.bg = RoundedRectangle(
            corner_radius=0.06,
            width=side,
            height=side,
            color=COLOR_BLACK,
            stroke_width=0,
            fill_color="#D3D3D3",
            fill_opacity=1.0,
        ).move_to(RIGHT * x_pos + UP * y_pos)

        hatch_group = VGroup()
        step = (2 * side) / (HATCH_LINES_SQUARE_GATE + 1)
        for i in range(1, HATCH_LINES_SQUARE_GATE + 1):
            c = -side + i * step
            x1 = max(-side / 2, c - side / 2)
            y1 = c - x1
            x2 = min(side / 2, c + side / 2)
            y2 = c - x2
            line = Line(
                RIGHT * (x_pos + x1) + UP * (y_pos + y1),
                RIGHT * (x_pos + x2) + UP * (y_pos + y2),
                stroke_width=2.0,
                color=color,
                stroke_opacity=1.0,
            )
            hatch_group.add(line)

        self.box = RoundedRectangle(
            corner_radius=0.06,
            width=side,
            height=side,
            color=COLOR_BLACK,
            stroke_width=WIRE_THICKNESS * side,
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
        self.dot = Dot(point=RIGHT * x_pos + UP * control_y, radius=(WIRE_THICKNESS * 2.0) / 100, color=WHITE)
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
        self.dot = Dot(point=RIGHT * x_pos + UP * control_y, radius=(WIRE_THICKNESS * 2.0) / 100, color=WHITE)
        min_y = min([t[0] for t in targets_info])
        self.cline = Line(self.dot.get_center(), RIGHT * x_pos + UP * min_y, stroke_width=WIRE_THICKNESS, color=WHITE).set_z_index(-1)
        self.add(self.dot, self.cline)

        self.gates = []
        for (y_pos, label_str, color, is_hatched) in targets_info:
            gate = HatchedGate(label_str, x_pos, y_pos, color) if is_hatched else QuantumGate(label_str, x_pos, y_pos, color)
            gate.set_z_index(1, family=False)
            self.gates.append(gate)
            self.add(gate)


def build_interleaved_wire_positions(n_main):
    wire_y_positions = []
    y_cursor = 3.0
    for i in range(n_main):
        wire_y_positions.append(y_cursor)
        y_cursor -= WIRE_SPACING_MAIN_COPY
        wire_y_positions.append(y_cursor)
        if i < n_main - 1:
            y_cursor -= WIRE_SPACING_PAIR

    y_top = max(wire_y_positions)
    y_bottom = min(wire_y_positions)
    y_center = 0.5 * (y_top + y_bottom)
    return [y - y_center for y in wire_y_positions]


class QFTCircuitPNG(Scene):
    def construct(self):
        set_text_defaults()

        n_main = 4
        qft_slots = 7
        wire_spacing = 1.35

        wire_edge_margin = MIN_GATE_GAP
        total_width = qft_slots * GATE_CENTER_SPACING + 2 * wire_edge_margin
        wire_start_x = CIRCUIT_CENTER_X - total_width / 2

        first_x = wire_start_x + wire_edge_margin + GATE_WIDTH / 2
        qft_x = [first_x + i * GATE_CENTER_SPACING for i in range(qft_slots)]
        y_positions = [2.1 - i * wire_spacing for i in range(n_main)]

        wires = []
        for i in range(n_main):
            wire = QuantumWire(
                y_pos=y_positions[i],
                length=total_width,
                start_x=wire_start_x,
                idx=i,
                is_copy=False,
                display_idx=i + 1,
                label_scale=0.62,
            )
            wires.append(wire)

        qft_items = []
        palette = [COLOR_BROWN_1, COLOR_BROWN_2, COLOR_BROWN_3]
        for i in range(n_main):
            qft_items.append(QuantumGate("H", qft_x[2 * i], y_positions[i], COLOR_H))
            if i < n_main - 1:
                targets = [
                    (
                        y_positions[tj],
                        f"P_{{{i + 1}{tj + 1}}}",
                        palette[min(tj - i - 1, 2)],
                        False,
                    )
                    for tj in range(i + 1, n_main)
                ]
                qft_items.append(PhaseTrain(qft_x[2 * i + 1], y_positions[i], targets))

        circuit = VGroup(*wires, *qft_items)
        circuit.scale(1.5)
        self.add(circuit)
        self.wait(1 / 30)


class DTCircuitPNG(Scene):
    def construct(self):
        set_text_defaults()

        n_main = 4
        n_total = 2 * n_main
        dt_slots, copy_slots = 7, 3
        section_gap = 1.0
        copy_spacing = GATE_CENTER_SPACING * 1.6
        wire_edge_margin = MIN_GATE_GAP

        total_width = dt_slots * GATE_CENTER_SPACING + section_gap + copy_slots * copy_spacing + 2 * wire_edge_margin
        wire_start_x = CIRCUIT_CENTER_X - total_width / 2

        first_x = wire_start_x + wire_edge_margin + GATE_WIDTH / 2
        dt_x = [first_x + i * GATE_CENTER_SPACING for i in range(dt_slots)]
        copy_x = [dt_x[-1] + GATE_CENTER_SPACING + section_gap + i * copy_spacing for i in range(copy_slots)]

        wire_y_positions = build_interleaved_wire_positions(n_main)

        def main_y(i):
            return wire_y_positions[2 * i]

        def copy_y(i):
            return wire_y_positions[2 * i + 1]

        wires = []
        for wi in range(n_total):
            wires.append(QuantumWire(wire_y_positions[wi], total_width, wire_start_x, wi, is_copy=(wi % 2 == 1)))

        dt_items = []
        palette = [COLOR_BROWN_1, COLOR_BROWN_2, COLOR_BROWN_3]
        for i in range(n_main):
            hd = HatchedGate("H_d", dt_x[2 * i], main_y(i), COLOR_HD)
            dt_items.append(hd)
            if i < n_main - 1:
                targets = [
                    (
                        main_y(tj),
                        f"R_{{{tj + 1}{i + 2}}}",
                        palette[min(i - tj, 2)],
                        True,
                    )
                    for tj in range(i + 1)
                ]
                dt_items.append(InvertedPhaseTrain(dt_x[2 * i + 1], main_y(i + 1), targets))

        copy_items = []
        for j in range(n_main - 1):
            targets = [
                (
                    main_y(ti),
                    f"R_{{{ti + 1}{j + 1}}}",
                    palette[min(ti - j - 1, 2)],
                    True,
                )
                for ti in range(j + 1, n_main)
            ]
            copy_items.append(PhaseTrain(copy_x[j], copy_y(j), targets))

        circuit = VGroup(*wires, *dt_items, *copy_items)
        circuit.scale(1.2)
        self.add(circuit)
        self.wait(1 / 30)


class ZTCircuitPNG(Scene):
    def construct(self):
        set_text_defaults()

        n_main = 4
        n_total = 2 * n_main
        dt_slots, copy_slots, qft_slots = 7, 3, 7
        section_gap = 0.5
        wire_edge_margin = MIN_GATE_GAP

        total_width = (dt_slots + copy_slots + qft_slots) * GATE_CENTER_SPACING + 2 * section_gap + 2 * wire_edge_margin
        wire_start_x = CIRCUIT_CENTER_X - total_width / 2

        first_x = wire_start_x + wire_edge_margin + GATE_WIDTH / 2
        dt_x = [first_x + i * GATE_CENTER_SPACING for i in range(dt_slots)]
        copy_x = [dt_x[-1] + GATE_CENTER_SPACING + section_gap + i * GATE_CENTER_SPACING for i in range(copy_slots)]
        qft_x = [copy_x[-1] + GATE_CENTER_SPACING + section_gap + i * GATE_CENTER_SPACING for i in range(qft_slots)]

        wire_y_positions = build_interleaved_wire_positions(n_main)

        def main_y(i):
            return wire_y_positions[2 * i]

        def copy_y(i):
            return wire_y_positions[2 * i + 1]

        wires = []
        for wi in range(n_total):
            wires.append(QuantumWire(wire_y_positions[wi], total_width, wire_start_x, wi, is_copy=(wi % 2 == 1)))

        palette = [COLOR_BROWN_1, COLOR_BROWN_2, COLOR_BROWN_3]

        dt_items = []
        for i in range(n_main):
            hd = HatchedGate("H_d", dt_x[2 * i], main_y(i), COLOR_HD)
            dt_items.append(hd)
            if i < n_main - 1:
                targets = [
                    (
                        main_y(tj),
                        f"R_{{{tj + 1}{i + 2}}}",
                        palette[min(i - tj, 2)],
                        True,
                    )
                    for tj in range(i + 1)
                ]
                dt_items.append(InvertedPhaseTrain(dt_x[2 * i + 1], main_y(i + 1), targets))

        copy_items = []
        for j in range(n_main - 1):
            targets = [
                (
                    main_y(ti),
                    f"R_{{{ti + 1}{j + 1}}}",
                    palette[min(ti - j - 1, 2)],
                    True,
                )
                for ti in range(j + 1, n_main)
            ]
            copy_items.append(PhaseTrain(copy_x[j], copy_y(j), targets))

        qft_items = []
        step = 0
        for orig_i in reversed(range(n_main)):
            mapped_i = (n_main - 1) - orig_i
            
            if orig_i < n_main - 1:
                targets = [
                    (
                        copy_y((n_main - 1) - orig_tj),
                        f"P_{{{orig_i + 1}{orig_tj + 1}}}",
                        palette[min(orig_tj - orig_i - 1, 2)],
                        False,
                    )
                    for orig_tj in range(orig_i + 1, n_main)
                ]
                qft_items.append(InvertedPhaseTrain(qft_x[step], copy_y(mapped_i), targets))
                step += 1
                
            qft_items.append(QuantumGate("H", qft_x[step], copy_y(mapped_i), COLOR_H))
            step += 1

        self.add(*wires, *dt_items, *copy_items, *qft_items)
        self.wait(1 / 30)
