from manim import *
import shutil
import os
from pathlib import Path

# Match the dark background used in tensorify.py
config.background_color = "#222831"

class SVDRMPSConversion(Scene):
    def construct(self):
        # Set default text colors to bright grey for contrast
        Text.set_default(color="#DDDDDD")
        MathTex.set_default(color="#DDDDDD")
        Tex.set_default(color="#DDDDDD")

        # Define theme colors
        COLOR_PURPLE = "#D24AD4"
        COLOR_YELLOW = "#E7C74A"
        COLOR_DARK_BLUE = "#0D1747"
        COLOR_LIGHT_BLUE = "#5DADE2"
        COLOR_BLACK = BLACK
        BOND_COLOR = COLOR_DARK_BLUE
        LEG_COLOR = BLACK
        U_WIDTH = 1.2
        SPLIT_GAP = 1.2 * U_WIDTH
        SPLIT_SHIFT = SPLIT_GAP / 2
        BASE_Y = 0.0
        S_SQUASH_Y = 0.06

        # Helpers
        def at_x(x):
            return RIGHT * x + UP * BASE_Y
        def make_tensor_rect(width, color, label_text=None):
            rect = RoundedRectangle(
                corner_radius=0.14,
                width=width,
                height=U_WIDTH,
                color=COLOR_BLACK,
                fill_opacity=0.85,
                fill_color=color,
                stroke_width=8,
            ).set_z_index(2)
            group = VGroup(rect)
            if label_text:
                label = MathTex(label_text, color=BLACK).set_z_index(3)
                group.add(label)
            return group

        def make_s_node(color, label_text="S"):
            square = Square(
                side_length=0.6,
                color=COLOR_BLACK,
                fill_opacity=0.9,
                fill_color=color,
                stroke_width=8,
            ).set_z_index(0)
            square.rotate(PI/4)
            group = VGroup(square)
            if label_text:
                label = MathTex(label_text, color=BLACK, font_size=32).set_z_index(0)
                group.add(label)
            return group

        def make_bond(start_point, end_point, split=True):
            thickness = 16 if split else 5
            bond_color = BOND_COLOR if split else COLOR_LIGHT_BLUE
            return Line(
                start=start_point,
                end=end_point,
                stroke_width=thickness,
                color=bond_color,
            ).set_z_index(1)

        def make_leg(anchor_point, label_text):
            line = Line(
                start=anchor_point,
                end=anchor_point + UP * 0.7,
                stroke_width=7,
                color=LEG_COLOR,
            ).set_z_index(1)
            tip = Circle(radius=0.07, color=LEG_COLOR, fill_opacity=1).move_to(line.get_end()).set_z_index(2)
            label = MathTex(label_text, color="#DDDDDD").next_to(tip, UP, buff=0.1).set_z_index(2)
            return VGroup(line, tip, label)

        def split_x_positions(parent_tensor, left_width, right_width):
            parent_body = parent_tensor[0]
            parent_left = parent_body.get_left()[0]
            left_temp_x = parent_left + left_width / 2
            right_temp_x = parent_left + left_width + right_width / 2
            left_final_x = left_temp_x - SPLIT_SHIFT
            right_final_x = right_temp_x + SPLIT_SHIFT
            return left_temp_x, right_temp_x, left_final_x, right_final_x

        # STEP 0: Initial full tensor T
        # T_body width = 4.8. Centered at 0. Range: [-2.4, 2.4]
        # Legs are spaced out evenly at x = -1.8, -0.6, 0.6, 1.8
        t_body = make_tensor_rect(4.8, COLOR_PURPLE, "T")
        t_body.move_to(at_x(0))
        
        leg_xs = [-1.8, -0.6, 0.6, 1.8]
        leg_anchor_y = t_body[0].get_top()[1]
        legs = VGroup(*[make_leg(RIGHT * x + UP * leg_anchor_y, f"n_{i+1}") for i, x in enumerate(leg_xs)])
        
        self.play(FadeIn(t_body), Create(legs))
        self.wait(1)

        leg1, leg2, leg3, leg4 = legs[0], legs[1], legs[2], legs[3]
        left_mobjects = []

        # ==================== STEP 1 ====================
        # u1 covers [-2.4, -1.2]. Center = -1.8. Width = 1.2.
        # v1 covers [-1.2, 2.4]. Center = 0.6. Width = 3.6.
        u1_temp_x, v1_temp_x, u1_x, v1_x = split_x_positions(t_body, U_WIDTH, 3.6)
        u1_temp = make_tensor_rect(U_WIDTH, COLOR_YELLOW, "U_1").move_to(at_x(u1_temp_x))
        v1_temp = make_tensor_rect(3.6, COLOR_PURPLE, "V_1").move_to(at_x(v1_temp_x))
        
        # Old tensor morphs into V1 while U1 is pulled out as a copy
        self.play(
            TransformFromCopy(t_body, u1_temp),
            Transform(t_body, v1_temp),
            run_time=0.5,
        )
        
        # Shift apart (shorter fixed split gap)
        u1 = make_tensor_rect(U_WIDTH, COLOR_YELLOW, "U_1").move_to(at_x(u1_x))
        v1 = make_tensor_rect(3.6, COLOR_PURPLE, "V_1").move_to(at_x(v1_x))
        s1_x = 0.5 * (u1[0].get_right()[0] + v1[0].get_left()[0])
        s1_target = make_s_node(COLOR_DARK_BLUE, "S_1").move_to(at_x(s1_x))
        s1 = s1_target.copy().stretch(S_SQUASH_Y, dim=1)
        self.add(s1)
        
        b1L = make_bond(u1[0].get_right(), s1[0].get_left(), split=True)
        b1R = make_bond(s1[0].get_right(), v1[0].get_left(), split=True)
        
        self.play(
            ReplacementTransform(u1_temp, u1),
            leg1.animate.shift(LEFT * SPLIT_SHIFT),
            Transform(t_body, v1),
            leg2.animate.shift(RIGHT * SPLIT_SHIFT),
            leg3.animate.shift(RIGHT * SPLIT_SHIFT),
            leg4.animate.shift(RIGHT * SPLIT_SHIFT),
            Transform(s1, s1_target), Create(b1L), Create(b1R),
            run_time=1.0,
            rate_func=smooth,
        )
        v1 = t_body
        self.wait(0.5)

        # Truncate S1
        s1_light = make_s_node(COLOR_LIGHT_BLUE, "S_1").move_to(s1)
        b1L_thin = make_bond(u1[0].get_right(), s1[0].get_left(), split=False)
        b1R_thin = make_bond(s1[0].get_right(), v1[0].get_left(), split=False)

        self.play(
            Transform(s1, s1_light),
            Transform(b1L, b1L_thin),
            Transform(b1R, b1R_thin),
            run_time=0.6
        )

        # Collapse S1 into a bond
        b1_full = make_bond(u1[0].get_right(), v1[0].get_left(), split=False)
        chi1 = MathTex(r"\chi_1", color=COLOR_LIGHT_BLUE).next_to(b1_full, UP, buff=0.15).set_z_index(3)
        s1_collapsed = s1.copy().stretch(S_SQUASH_Y, dim=1).move_to(b1_full.get_center())

        self.play(
            ReplacementTransform(b1L, b1_full),
            FadeOut(b1R),
            Transform(s1, s1_collapsed),
            FadeIn(chi1),
            run_time=1
        )
        self.remove(s1)
        self.wait(0.5)
        
        left_mobjects.extend([u1, leg1, b1_full, chi1])

        # ==================== STEP 2 ====================
        # v1 currently covers [0.0, 3.6]. Center = 1.8. Width = 3.6.
        # Splitting v1 into:
        # u2 covers [0.0, 1.2]. Center = 0.6. Width 1.2.
        # v2 covers [1.2, 3.6]. Center = 2.4. Width 2.4.
        u2_temp_x, v2_temp_x, u2_x, v2_x = split_x_positions(v1, U_WIDTH, 2.4)
        u2_temp = make_tensor_rect(U_WIDTH, COLOR_YELLOW, "U_2").move_to(at_x(u2_temp_x))
        v2_temp = make_tensor_rect(2.4, COLOR_PURPLE, "V_2").move_to(at_x(v2_temp_x))
        
        self.play(
            TransformFromCopy(v1, u2_temp),
            Transform(v1, v2_temp),
            run_time=0.5,
        )
        
        # Shift apart (shorter fixed split gap)
        u2 = make_tensor_rect(U_WIDTH, COLOR_YELLOW, "U_2").move_to(at_x(u2_x))
        v2 = make_tensor_rect(2.4, COLOR_PURPLE, "V_2").move_to(at_x(v2_x))
        s2_x = 0.5 * (u2[0].get_right()[0] + v2[0].get_left()[0])
        s2_target = make_s_node(COLOR_DARK_BLUE, "S_2").move_to(at_x(s2_x))
        s2 = s2_target.copy().stretch(S_SQUASH_Y, dim=1)
        self.add(s2)
        
        b2L = make_bond(u2[0].get_right(), s2[0].get_left(), split=True)
        b2R = make_bond(s2[0].get_right(), v2[0].get_left(), split=True)
        
        shift_Left_Step2 = [mob.animate.shift(LEFT * SPLIT_SHIFT) for mob in left_mobjects]
        
        self.play(
            ReplacementTransform(u2_temp, u2),
            leg2.animate.shift(LEFT * SPLIT_SHIFT),
            Transform(v1, v2),
            leg3.animate.shift(RIGHT * SPLIT_SHIFT),
            leg4.animate.shift(RIGHT * SPLIT_SHIFT),
            Transform(s2, s2_target), Create(b2L), Create(b2R),
            *shift_Left_Step2,
            run_time=1.0,
            rate_func=smooth,
        )
        v2 = v1
        self.wait(0.5)

        # Truncate S2
        s2_light = make_s_node(COLOR_LIGHT_BLUE, "S_2").move_to(s2)
        b2L_thin = make_bond(u2[0].get_right(), s2[0].get_left(), split=False)
        b2R_thin = make_bond(s2[0].get_right(), v2[0].get_left(), split=False)

        self.play(
            Transform(s2, s2_light),
            Transform(b2L, b2L_thin),
            Transform(b2R, b2R_thin),
            run_time=0.6
        )

        # Collapse S2 into a bond
        b2_full = make_bond(u2[0].get_right(), v2[0].get_left(), split=False)
        chi2 = MathTex(r"\chi_2", color=COLOR_LIGHT_BLUE).next_to(b2_full, UP, buff=0.15).set_z_index(3)
        s2_collapsed = s2.copy().stretch(S_SQUASH_Y, dim=1).move_to(b2_full.get_center())

        self.play(
            ReplacementTransform(b2L, b2_full),
            FadeOut(b2R),
            Transform(s2, s2_collapsed),
            FadeIn(chi2),
            run_time=1
        )
        self.remove(s2)
        self.wait(0.5)
        
        left_mobjects.extend([u2, leg2, b2_full, chi2])

        # ==================== STEP 3 ====================
        # v2 currently covers [2.4, 4.8]. Center = 3.6. Width = 2.4.
        # Splitting v2 into:
        # u3 covers [2.4, 3.6]. Center = 3.0. Width 1.2.
        # v3 covers [3.6, 4.8]. Center = 4.2. Width 1.2.
        u3_temp_x, v3_temp_x, u3_x, v3_x = split_x_positions(v2, U_WIDTH, U_WIDTH)
        u3_temp = make_tensor_rect(U_WIDTH, COLOR_YELLOW, "U_3").move_to(at_x(u3_temp_x))
        v3_temp = make_tensor_rect(U_WIDTH, COLOR_PURPLE, "V_3").move_to(at_x(v3_temp_x))
        
        self.play(
            TransformFromCopy(v2, u3_temp),
            Transform(v2, v3_temp),
            run_time=0.5,
        )
        
        # Shift apart (shorter fixed split gap)
        u3 = make_tensor_rect(U_WIDTH, COLOR_YELLOW, "U_3").move_to(at_x(u3_x))
        v3 = make_tensor_rect(U_WIDTH, COLOR_PURPLE, "V_3").move_to(at_x(v3_x))
        s3_x = 0.5 * (u3[0].get_right()[0] + v3[0].get_left()[0])
        s3_target = make_s_node(COLOR_DARK_BLUE, "S_3").move_to(at_x(s3_x))
        s3 = s3_target.copy().stretch(S_SQUASH_Y, dim=1)
        self.add(s3)
        
        b3L = make_bond(u3[0].get_right(), s3[0].get_left(), split=True)
        b3R = make_bond(s3[0].get_right(), v3[0].get_left(), split=True)
        
        shift_Left_Step3 = [mob.animate.shift(LEFT * SPLIT_SHIFT) for mob in left_mobjects]
        
        self.play(
            ReplacementTransform(u3_temp, u3),
            leg3.animate.shift(LEFT * SPLIT_SHIFT),
            Transform(v2, v3),
            leg4.animate.shift(RIGHT * SPLIT_SHIFT),
            Transform(s3, s3_target), Create(b3L), Create(b3R),
            *shift_Left_Step3,
            run_time=1.0,
            rate_func=smooth,
        )
        v3 = v2
        self.wait(0.5)

        # Truncate S3
        s3_light = make_s_node(COLOR_LIGHT_BLUE, "S_3").move_to(s3)
        b3L_thin = make_bond(u3[0].get_right(), s3[0].get_left(), split=False)
        b3R_thin = make_bond(s3[0].get_right(), v3[0].get_left(), split=False)

        self.play(
            Transform(s3, s3_light),
            Transform(b3L, b3L_thin),
            Transform(b3R, b3R_thin),
            run_time=0.6
        )

        # Collapse S3 into a bond
        b3_full = make_bond(u3[0].get_right(), v3[0].get_left(), split=False)
        chi3 = MathTex(r"\chi_3", color=COLOR_LIGHT_BLUE).next_to(b3_full, UP, buff=0.15).set_z_index(3)
        s3_collapsed = s3.copy().stretch(S_SQUASH_Y, dim=1).move_to(b3_full.get_center())

        self.play(
            ReplacementTransform(b3L, b3_full),
            FadeOut(b3R),
            Transform(s3, s3_collapsed),
            FadeIn(chi3),
            run_time=1
        )
        self.remove(s3)
        self.wait(1.5)

        # Auto-copy to assets folder as in tensorify.py
        try:
            fw = getattr(self.renderer, "file_writer", None)
            movie_path = None
            if fw is not None:
                movie_path = getattr(fw, "movie_file_path", None) or getattr(fw, "movie_path", None)

            assets_dir = Path(__file__).parent / "assets"
            assets_dir.mkdir(parents=True, exist_ok=True)

            if movie_path and os.path.exists(movie_path):
                dest = assets_dir / Path(movie_path).name
                shutil.copy(movie_path, dest)
            else:
                search_root = Path(__file__).parent / "media" / "videos"
                candidates = list(search_root.rglob("*SVDRMPSConversion*.mp4"))
                if candidates:
                    latest = max(candidates, key=lambda p: p.stat().st_mtime)
                    dest = assets_dir / latest.name
                    shutil.copy(latest, dest)
        except Exception as e:
            print("Post-render copy to animations/assets failed:", e)


