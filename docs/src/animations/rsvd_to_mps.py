from manim import *
import shutil
import os
from pathlib import Path

# Match the dark background used in tensorify.py
config.background_color = "#222831"

class RSVDRMPSConversion(Scene):
    def construct(self):
        # Set default text colors to bright grey for contrast
        Text.set_default(color="#DDDDDD")
        MathTex.set_default(color="#DDDDDD")
        Tex.set_default(color="#DDDDDD")

        # Define theme colors
        COLOR_PURPLE = "#D24AD4"
        COLOR_YELLOW = "#E7C74A"
        COLOR_DARK_BLUE = "#0D1747"
        COLOR_LIGHT_BLUE = "#297DB5"
        COLOR_BLACK = BLACK
        BOND_COLOR = COLOR_LIGHT_BLUE
        LEG_COLOR = BLACK
        
        # Site dimensions and spacing
        U_WIDTH = 1.2
        SPLIT_GAP = 1.2  # Defined by user: same as site width
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
            ).set_z_index(3)
            group = VGroup(rect)
            if label_text:
                label = MathTex(label_text, color=BLACK).set_z_index(4)
                group.add(label)
            return group

        def make_s_node(color, label_text="S"):
            square = Square(
                side_length=0.6,
                color=COLOR_BLACK,
                fill_opacity=0.9,
                fill_color=color,
                stroke_width=8,
            ).set_z_index(2)
            square.rotate(PI/4)
            group = VGroup(square)
            if label_text:
                label = MathTex(label_text, color=BLACK, font_size=32).set_z_index(3)
                group.add(label)
            return group

        def make_bond(start_point, end_point, split=True):
            thickness = 24 if split else 5
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
            label = MathTex(label_text, color="#DDDDDD").next_to(tip, UP, buff=0.1).set_z_index(3)
            return VGroup(line, tip, label)

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

        # ==================== STEP 1 ====================
        # Split T into U_1 (Yellow, Left half) and V_1 (Purple, Right half)
        # Shift slightly so their gap is 1.2
        # U_1 center = -1.8, V_1 center = 1.8
        
        u_temp = make_tensor_rect(2.4, COLOR_YELLOW, "U_1").move_to(at_x(-1.2)) # Just left of center
        v_temp = make_tensor_rect(2.4, COLOR_PURPLE, "V_1").move_to(at_x(1.2)) # Just right of center

        self.play(
            TransformFromCopy(t_body, u_temp),
            Transform(t_body, v_temp),
            run_time=0.5,
        )

        u_final = make_tensor_rect(2.4, COLOR_YELLOW, "U_1").move_to(at_x(-1.8))
        v_final = make_tensor_rect(2.4, COLOR_PURPLE, "V_1").move_to(at_x(1.8))

        s_mid_target = make_s_node(COLOR_LIGHT_BLUE, "S_2").move_to(at_x(0))
        s_mid = s_mid_target.copy().stretch(S_SQUASH_Y, dim=1)
        self.add(s_mid)

        b_mid = make_bond(u_final[0].get_right(), v_final[0].get_left(), split=True)

        self.play(
            ReplacementTransform(u_temp, u_final),
            leg1.animate.shift(LEFT * 0.6), # Move from -1.8 to -2.4 
            leg2.animate.shift(LEFT * 0.6), # Move from -0.6 to -1.2
            Transform(t_body, v_final),
            leg3.animate.shift(RIGHT * 0.6), # Move from 0.6 to 1.2
            leg4.animate.shift(RIGHT * 0.6), # Move from 1.8 to 2.4
            Transform(s_mid, s_mid_target),
            Create(b_mid),
            run_time=1.0,
            rate_func=smooth,
        )
        self.wait(0.5)

        # Truncate S_mid
        s_mid_light = make_s_node(COLOR_LIGHT_BLUE, "S_2").move_to(s_mid)
        b_mid_thin = make_bond(u_final[0].get_right(), v_final[0].get_left(), split=False)

        self.play(
            Transform(s_mid, s_mid_light),
            Transform(b_mid, b_mid_thin),
            run_time=0.6
        )

        # Collapse S_mid into bond chi_2
        chi_mid = MathTex(r"\chi_2", color=COLOR_LIGHT_BLUE).next_to(b_mid_thin, UP, buff=0.15).set_z_index(5)
        s_mid_collapsed = s_mid.copy().stretch(S_SQUASH_Y, dim=1).move_to(b_mid_thin.get_center())

        self.play(
            Transform(s_mid, s_mid_collapsed),
            FadeIn(chi_mid),
            run_time=0.5
        )
        self.remove(s_mid)
        self.wait(0.5)

        # ==================== STEP 2 ====================
        # Recursive splitting into U2, V2, U3, V3
        # U_1 (-1.8) -> U_2 (-2.4, Yellow) and V_2 (-1.2, Purple)
        # V_1 (1.8)  -> U_3 (1.2, Yellow) and V_3 (2.4, Purple)
        u2_temp = make_tensor_rect(1.2, COLOR_YELLOW, "U_2").move_to(at_x(-2.4))
        v2_temp = make_tensor_rect(1.2, COLOR_PURPLE, "V_2").move_to(at_x(-1.2))
        
        u3_temp = make_tensor_rect(1.2, COLOR_YELLOW, "U_3").move_to(at_x(1.2))
        v3_temp = make_tensor_rect(1.2, COLOR_PURPLE, "V_3").move_to(at_x(2.4))

        self.play(
            TransformFromCopy(u_final, u2_temp),
            Transform(u_final, v2_temp),
            TransformFromCopy(t_body, v3_temp),
            Transform(t_body, u3_temp),
            run_time=0.5,
        )

        # Move to spread out to final positions
        # Gap needs to be exactly 1.2 everywhere
        u2_final = make_tensor_rect(1.2, COLOR_YELLOW, "U_2").move_to(at_x(-3.6))
        v2_final = make_tensor_rect(1.2, COLOR_PURPLE, "V_2").move_to(at_x(-1.2))
        u3_final = make_tensor_rect(1.2, COLOR_YELLOW, "U_3").move_to(at_x(1.2))
        v3_final = make_tensor_rect(1.2, COLOR_PURPLE, "V_3").move_to(at_x(3.6))

        s_L_target = make_s_node(COLOR_LIGHT_BLUE, "S_1").move_to(at_x(-2.4))
        s_L = s_L_target.copy().stretch(S_SQUASH_Y, dim=1)
        self.add(s_L)

        s_R_target = make_s_node(COLOR_LIGHT_BLUE, "S_3").move_to(at_x(2.4))
        s_R = s_R_target.copy().stretch(S_SQUASH_Y, dim=1)
        self.add(s_R)

        b_L = make_bond(u2_final[0].get_right(), v2_final[0].get_left(), split=True)
        b_R = make_bond(u3_final[0].get_right(), v3_final[0].get_left(), split=True)

        # Updaters for the middle bond since its anchors change
        def update_mid_bond(m):
            m.become(make_bond(v2_final[0].get_right(), u3_final[0].get_left(), split=False).set_z_index(1))
        def update_chi_mid(m):
            m.next_to(b_mid_thin, UP, buff=0.15).set_z_index(5)
            
        b_mid_thin.add_updater(update_mid_bond)
        chi_mid.add_updater(update_chi_mid)

        self.play(
            ReplacementTransform(u2_temp, u2_final),
            Transform(u_final, v2_final), 
            Transform(t_body, u3_final),
            ReplacementTransform(v3_temp, v3_final),
            
            leg1.animate.shift(LEFT * 1.2), # -2.4 -> -3.6
            # leg2 is at -1.2, stays at -1.2 (moves relative to U_1, but V_2 center is -1.2)
            # leg3 is at 1.2, stays at 1.2
            leg4.animate.shift(RIGHT * 1.2), # 2.4 -> 3.6

            Transform(s_L, s_L_target),
            Transform(s_R, s_R_target),
            Create(b_L), Create(b_R),
            run_time=1.0,
            rate_func=smooth,
        )
        b_mid_thin.clear_updaters()
        chi_mid.clear_updaters()
        
        self.wait(0.5)

        # Truncate S_L and S_R
        s_L_light = make_s_node(COLOR_LIGHT_BLUE, "S_1").move_to(s_L)
        s_R_light = make_s_node(COLOR_LIGHT_BLUE, "S_3").move_to(s_R)

        b_L_thin = make_bond(u2_final[0].get_right(), v2_final[0].get_left(), split=False)
        b_R_thin = make_bond(u3_final[0].get_right(), v3_final[0].get_left(), split=False)

        self.play(
            Transform(s_L, s_L_light),
            Transform(s_R, s_R_light),
            Transform(b_L, b_L_thin),
            Transform(b_R, b_R_thin),
            run_time=0.6
        )

        # Collapse S_L and S_R into bonds
        chi_L = MathTex(r"\chi_1", color=COLOR_LIGHT_BLUE).next_to(b_L_thin, UP, buff=0.15).set_z_index(5)
        s_L_collapsed = s_L.copy().stretch(S_SQUASH_Y, dim=1).move_to(b_L_thin.get_center())

        chi_R = MathTex(r"\chi_3", color=COLOR_LIGHT_BLUE).next_to(b_R_thin, UP, buff=0.15).set_z_index(5)
        s_R_collapsed = s_R.copy().stretch(S_SQUASH_Y, dim=1).move_to(b_R_thin.get_center())

        self.play(
            Transform(s_L, s_L_collapsed),
            Transform(s_R, s_R_collapsed),
            FadeIn(chi_L), FadeIn(chi_R),
            run_time=0.5
        )
        self.remove(s_L, s_R)
        self.wait(1.5)

        # Auto-copy to assets folder
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
                candidates = list(search_root.rglob("*RSVDRMPSConversion*.mp4"))
                if candidates:
                    latest = max(candidates, key=lambda p: p.stat().st_mtime)
                    dest = assets_dir / latest.name
                    shutil.copy(latest, dest)
        except Exception as e:
            print("Post-render copy to animations/assets failed:", e)
