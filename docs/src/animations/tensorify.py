from manim import *
import shutil
import os
from pathlib import Path

# Set the background color to the requested dark shade.
config.background_color = "#222831"

class TensorifyVector(ThreeDScene):
    def construct(self):
        # Use a light text color for readability on dark background.
        Text.set_default(color="#DDDDDD")
        MathTex.set_default(color="#DDDDDD")
        Tex.set_default(color="#DDDDDD")

        # 1. Start with a 1D array of 8 numbers
        
        # Define 8 elements
        vector_elements = VGroup(*[MathTex(f"x_{i}") for i in range(8)])
        vector_elements.arrange(RIGHT, buff=0.7)
        
        # Draw a rounded box around each element to make it look like an array
        boxes = VGroup(*[RoundedRectangle(corner_radius=0.15, width=1.0, height=0.8, color=BLUE, stroke_width=3) for _ in range(8)])
        for box, elem in zip(boxes, vector_elements):
            box.move_to(elem.get_center())
            
        vector_group = VGroup()
        for i in range(8):
            cell = VGroup(boxes[i], vector_elements[i])
            vector_group.add(cell)
            
        # We don't fix vector_group in frame yet because we want to move them in 3D space
        self.play(Write(vector_group), run_time=1)
        
        # 2. Transform the 1D array into a 3D tensor (2x2x2)
        # Set up 3D camera
        self.move_camera(phi=65 * DEGREES, theta=-55 * DEGREES, run_time=1.5)
        
        # We need to calculate the 2x2x2 grid positions
        grid_positions = []
        for k in range(2):
            for j in range(2):
                for i in range(2):
                    # Shift everything LEFT by 1.5 and DOWN by 1.0 to center the 3D object lower
                    pos = LEFT * 1.5 + DOWN * 1.0 + RIGHT * (i - 0.5) * 1.5 + DOWN * (j - 0.5) * 1.5 + OUT * (k - 0.5) * 1.5
                    grid_positions.append(pos)
        
        # Create axes arrows for the 3 indices (drawn from corner), shifted left
        axes_origin = LEFT * 4.5 + DOWN * 3.
        axis_i = Arrow(start=axes_origin, end=axes_origin + RIGHT * 1.5, color=RED, stroke_width=4)
        axis_j = Arrow(start=axes_origin, end=axes_origin + DOWN * 1.5, color=GREEN, stroke_width=4)
        axis_k = Arrow(start=axes_origin, end=axes_origin + OUT * 1.5, color=ORANGE, stroke_width=4)
        
        # Add index labels to the axes
        label_i = MathTex("i_3").next_to(axis_i.get_end(), RIGHT, buff=0.1).set_color(RED)
        label_j = MathTex("i_2").next_to(axis_j.get_end(), DOWN, buff=0.1).set_color(GREEN)
        label_k = MathTex("i_1").next_to(axis_k.get_end(), OUT, buff=0.1).set_color(ORANGE)
        
        axes_group = VGroup(axis_i, axis_j, axis_k, label_i, label_j, label_k)
        
        # The final tensor must be a rounded rectangle with all legs pointing up
        # It should appear on the right side, with a thick black border.
        final_tensor = RoundedRectangle(corner_radius=0.4, width=3.6, height=1.5, color=BLACK, fill_opacity=0.4, fill_color="#D24AD4", stroke_width=8)
        final_tensor.move_to(RIGHT * 3.5 + DOWN * 0.5)
        # Put tensor label on top
        tensor_text = Text("T").move_to(final_tensor.get_center())
        final_tensor.set_z_index(5)
        tensor_text.set_z_index(6)
        
        # 3 Shorter legs pointing UPwards as requested
        # Sending legs "behind" the tensor by making them slightly shorter and keeping standard z-order
        leg_i1 = Line(final_tensor.get_top() + LEFT * 1.0, final_tensor.get_top() + LEFT * 1.0 + UP * 0.8, stroke_width=8, color=ORANGE)
        leg_i2 = Line(final_tensor.get_top(), final_tensor.get_top() + UP * 0.8, stroke_width=8, color=GREEN)
        leg_i3 = Line(final_tensor.get_top() + RIGHT * 1.0, final_tensor.get_top() + RIGHT * 1.0 + UP * 0.8, stroke_width=8, color=RED)
        
        # Add small circles at the top tip of the legs
        circ_i1 = Circle(radius=0.1, color=ORANGE, fill_opacity=1).move_to(leg_i1.get_end())
        circ_i2 = Circle(radius=0.1, color=GREEN, fill_opacity=1).move_to(leg_i2.get_end())
        circ_i3 = Circle(radius=0.1, color=RED, fill_opacity=1).move_to(leg_i3.get_end())
        
        # Labels for the legs (user wants indices i1, i2, i3 appearing on the right of the legs/object)
        leg_label_i1 = MathTex("i_1").next_to(circ_i1, RIGHT, buff=0.1).set_color(ORANGE)
        leg_label_i2 = MathTex("i_2").next_to(circ_i2, RIGHT, buff=0.1).set_color(GREEN)
        leg_label_i3 = MathTex("i_3").next_to(circ_i3, RIGHT, buff=0.1).set_color(RED)
        leg_label_i1.set_z_index(2)
        leg_label_i2.set_z_index(2)
        leg_label_i3.set_z_index(2)
        
        # Separate the body and legs for the sequential animation
        tensor_body = VGroup(final_tensor, tensor_text)
        tensor_legs = VGroup(leg_i1, leg_i2, leg_i3, circ_i1, circ_i2, circ_i3, leg_label_i1, leg_label_i2, leg_label_i3)
        
        animations_p1 = []
        for idx in range(8):
            animations_p1.append(vector_group[idx].animate.move_to(grid_positions[idx]))
        
        animations_p1.append(FadeIn(axes_group))
        
        # We prepare the 2D overlay animations
        # Register them as fixed-in-frame, but keep them out of the scene until delayed reveal.
        self.add_fixed_in_frame_mobjects(tensor_legs, tensor_body)
        self.remove(tensor_legs, tensor_body)

        # Start cube rearrangement immediately; duration is controlled here.
        move_cube = AnimationGroup(*animations_p1, run_time=1.5)
        show_tensor = Succession(
            Wait(0.8),
            AnimationGroup(
                FadeIn(tensor_body),
                Create(tensor_legs),
                run_time=1.0,
            ),
        )

        self.play(move_cube, show_tensor)

        # Move tensor up to make room for the index-value narration below.
        self.play(
            tensor_body.animate.shift(UP * 1.2),
            tensor_legs.animate.shift(UP * 1.2),
            run_time=0.8,
        )

        def tensor_symbolic_tex():
            symbolic = MathTex(
                "T", "[", "i_1", ",", "i_2", ",", "i_3", "]", "=", "x_j"
            )
            symbolic[0].set_color(BLUE_A)
            symbolic[2].set_color(ORANGE)
            symbolic[4].set_color(GREEN)
            symbolic[6].set_color(RED)
            symbolic.next_to(tensor_body, DOWN, buff=1.35)
            return symbolic

        def tensor_assignment_tex(index_value, i1, i2, i3):
            assignment = MathTex(
                "T", "[", str(i1), ",", str(i2), ",", str(i3), "]", "=", f"x_{{{index_value}}}"
            )
            assignment[0].set_color(BLUE_A)
            assignment[2].set_color(ORANGE)
            assignment[4].set_color(GREEN)
            assignment[6].set_color(RED)
            assignment.next_to(tensor_body, DOWN, buff=1.35)
            return assignment

        # Start with symbolic relation and no highlighted cube entry.
        highlight_lift = 0.3
        equation_tex = tensor_symbolic_tex()
        self.add_fixed_in_frame_mobjects(equation_tex)
        self.play(Write(equation_tex), run_time=0.95)
        self.wait(0.5)

        previous_index = None
        for idx in range(8):
            i1 = idx // 4 + 1
            i2 = (idx % 4) // 2 + 1
            i3 = idx % 2 + 1

            new_equation_tex = tensor_assignment_tex(idx, i1, i2, i3)
            new_equation_tex.move_to(equation_tex.get_center())

            reset_previous = []
            if previous_index is not None:
                reset_previous = [
                    vector_elements[previous_index].animate.shift(IN * highlight_lift).set_color("#DDDDDD"),
                    boxes[previous_index].animate.shift(IN * highlight_lift).set_color(BLUE),
                ]

            self.play(
                *reset_previous,
                Transform(equation_tex, new_equation_tex),
                vector_elements[idx].animate.shift(OUT * highlight_lift).set_color(RED),
                boxes[idx].animate.shift(OUT * highlight_lift).set_color(RED),
                run_time=0.95,
            )

            previous_index = idx

        # End with symbolic relation and no highlighted cube entry.
        final_symbolic = tensor_symbolic_tex()
        final_symbolic.move_to(equation_tex.get_center())
        self.play(
            vector_elements[previous_index].animate.shift(IN * highlight_lift).set_color("#DDDDDD"),
            boxes[previous_index].animate.shift(IN * highlight_lift).set_color(BLUE),
            Transform(equation_tex, final_symbolic),
            run_time=0.95,
        )

        self.wait(1)

        # After rendering, attempt to copy the generated mp4 into the animations/assets folder
        try:
            # Try to get movie path from manim's file writer (works on most manim versions)
            fw = getattr(self.renderer, "file_writer", None)
            movie_path = None
            if fw is not None:
                movie_path = getattr(fw, "movie_file_path", None) or getattr(fw, "movie_path", None)

            # If we have a valid movie path, copy it. Otherwise, fall back to searching common output folders.
            assets_dir = Path(__file__).parent / "assets"
            assets_dir.mkdir(parents=True, exist_ok=True)

            if movie_path and os.path.exists(movie_path):
                dest = assets_dir / Path(movie_path).name
                shutil.copy(movie_path, dest)
            else:
                # Fallback: search for the newest mp4 produced by this scene
                search_root = Path(__file__).parent / "media" / "videos"
                candidates = list(search_root.rglob("*TensorifyVector*.mp4"))
                if candidates:
                    latest = max(candidates, key=lambda p: p.stat().st_mtime)
                    dest = assets_dir / latest.name
                    shutil.copy(latest, dest)
        except Exception as e:
            # Do not crash the scene if copy fails; print error for debugging
            print("Post-render copy to animations/assets failed:", e)
