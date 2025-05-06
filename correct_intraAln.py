#!/usr/bin/env python3
import os
import tkinter as tk
from PIL import Image  # Avoid using ImageTk
import io
import base64
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--img-dir", default="intraAln/preview", help="Directory containing the images")
parser.add_argument("--ref-img", default="doublet_8nm.png", help="Reference image path")
parser.add_argument("--polarity", default="filamentListPolarity.csv", help="CSV file for polarity mapping")
args = parser.parse_args()

IMAGE_DIR = args.img_dir  # Change this to your directory
REFERENCE_IMAGE_PATH = args.ref_img  # Change this to your reference image
POLARITY_FILE = args.polarity  # Change this as needed

# New color variables
COLOR_0 = "#0000FF"  # Blue
COLOR_1 = "#FF0000"  # Red

# Read polarity mapping from CSV file
polarity_map = {}
try:
    with open(POLARITY_FILE, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split(",")
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = parts[1].strip()
                    polarity_map[key] = value
except Exception as e:
    print(f"Error reading polarity file: {e}")

# Get list of PNG images from the directory (excluding the reference image)
images = [
    os.path.join(IMAGE_DIR, f)
    for f in os.listdir(IMAGE_DIR)
    if f.lower().endswith("aln.png") and f != os.path.basename(REFERENCE_IMAGE_PATH) and not f.lower().endswith("manual_aln.png")
]

images.sort()


class PhotoBrowser:
    def __init__(self, master):
        self.master = master
        self.scale_factor = 1.0
        self.manual_editors = []  # open manual correction windows
        self.ref_alpha = 102  # default ref image opacity (0-255)
        self.ref_black = 0    # default black cut off
        self.ref_white = 255  # default white cut off

        # Object Viewer attributes (independent image sizing)
        self.object_viewer_window = None
        self.ov_scale_factor = 1.0  # independent scale factor for Object Viewer
        self.ov_current_family = None  # currently displayed family

        # Create main window menu
        menubar = tk.Menu(master)
        view_menu = tk.Menu(menubar, tearoff=0)
        view_menu.add_command(label="Increase Image Size (+)", command=self.increase_image_size)
        view_menu.add_command(label="Decrease Image Size (-)", command=self.decrease_image_size)
        view_menu.add_command(label="Default Image Size (0)", command=self.default_image_size)
        view_menu.add_command(label="Object Viewer (v)", command=self.object_viewer)
        menubar.add_cascade(label="View", menu=view_menu)
        go_menu = tk.Menu(menubar, tearoff=0)
        go_menu.add_command(label="Select Image...", command=self.open_go_window)
        menubar.add_cascade(label="Go", menu=go_menu)
        master.config(menu=menubar)

        # Bind keys in the main window
        self.master.bind("<KeyPress-equal>", lambda event: self.increase_image_size())
        self.master.bind("<KeyPress-KP_Add>", lambda event: self.increase_image_size())
        self.master.bind("<KeyPress-minus>", lambda event: self.decrease_image_size())
        self.master.bind("<KeyPress-0>", lambda event: self.default_image_size())
        self.master.bind("<Left>", self.prev_image)
        self.master.bind("<Right>", self.next_image)
        self.master.bind("v", lambda event: self.object_viewer())

        self.current_index = 0

        self.left_frame = tk.Frame(master)
        self.left_frame.pack(side="left", padx=10, pady=10)
        self.right_frame = tk.Frame(master)
        self.right_frame.pack(side="left", padx=10, pady=10)
        self.image_frame = tk.Frame(self.left_frame)
        self.image_frame.pack()

        self.ref_photo = self.load_and_scale_image(REFERENCE_IMAGE_PATH)
        self.ref_label = tk.Label(self.image_frame, image=self.ref_photo)
        self.ref_label.pack(side="left", padx=5, pady=5)
        self.img_label = tk.Label(self.image_frame)
        self.img_label.pack(side="left", padx=5, pady=5)

        self.manual_button = tk.Button(self.left_frame, text="Manual Correct", command=self.manual_correct)
        self.manual_button.pack(pady=10)
        self.delete_manual_button = tk.Button(self.left_frame, text="Delete Manual", command=self.delete_manual)
        self.delete_manual_button.pack(pady=5)
        
        # Two separate labels for filename and polarity
        self.file_info_label = tk.Label(self.right_frame, text="", justify="left", font=("Arial", 12))
        self.file_info_label.pack(anchor="nw")
        self.polarity_label = tk.Label(self.right_frame, text="", justify="left", font=("Arial", 12))
        self.polarity_label.pack(anchor="nw")
        
        self.transform_frame = tk.Frame(self.right_frame)
        self.transform_frame.pack(anchor="nw", pady=10)

        self.show_image()

    def load_and_scale_image(self, path):
        try:
            pil_image = Image.open(path).convert("RGBA")
            if self.scale_factor != 1.0:
                new_size = (int(pil_image.width * self.scale_factor), int(pil_image.height * self.scale_factor))
                pil_image = pil_image.resize(new_size, Image.ANTIALIAS)
            buffer = io.BytesIO()
            pil_image.save(buffer, format="PNG")
            encoded = base64.b64encode(buffer.getvalue())
            return tk.PhotoImage(data=encoded.decode("ascii"))
        except Exception as e:
            print(f"Error loading and scaling image {path}: {e}")
            return None

    def load_and_scale_image_OV(self, path):
        try:
            pil_image = Image.open(path).convert("RGBA")
            if self.ov_scale_factor != 1.0:
                new_size = (int(pil_image.width * self.ov_scale_factor), int(pil_image.height * self.ov_scale_factor))
                pil_image = pil_image.resize(new_size, Image.ANTIALIAS)
            buffer = io.BytesIO()
            pil_image.save(buffer, format="PNG")
            encoded = base64.b64encode(buffer.getvalue())
            return tk.PhotoImage(data=encoded.decode("ascii"))
        except Exception as e:
            print(f"Error loading and scaling image (OV) {path}: {e}")
            return None

    def extract_key(self, filename):
        if filename.endswith("_aln.png"):
            return filename[:-8]
        return filename[:-4] if filename.endswith(".png") else filename

    def get_tbl_data(self, key):
        tbl_path = os.path.join("particles", key, "xform.tbl")
        try:
            with open(tbl_path, "r") as tbl_file:
                line = tbl_file.readline().strip()
                if line:
                    parts = line.split()
                    if len(parts) >= 6:
                        return parts[:6]
        except Exception:
            pass
        return None

    def get_manual_tbl_data(self, key):
        tbl_path = os.path.join("particles", key, "manual_xform.tbl")
        try:
            with open(tbl_path, "r") as tbl_file:
                line = tbl_file.readline().strip()
                if line:
                    parts = line.split()
                    if len(parts) >= 6:
                        return parts[:6]
        except Exception:
            pass
        return None

    def show_image(self):
        if images:
            image_path = images[self.current_index]
            file_name = os.path.basename(image_path)
            self.master.title(file_name)
            self.img_photo = self.load_and_scale_image(image_path)
            self.img_label.config(image=self.img_photo)
            self.ref_photo = self.load_and_scale_image(REFERENCE_IMAGE_PATH)
            self.ref_label.config(image=self.ref_photo)
            key = self.extract_key(file_name)
            polarity = polarity_map.get(key, "N/A")
            fg_color = "black"
            if polarity == "1":
                fg_color = COLOR_1
            elif polarity == "0":
                fg_color = COLOR_0
            self.file_info_label.config(text=file_name, fg="black")
            self.polarity_label.config(text=f"Polarity: {polarity}", fg=fg_color)
            for widget in self.transform_frame.winfo_children():
                widget.destroy()
            labels = ["x-shift", "y-shift", "z-shift", "z-rot", "x-rot", "z-rot-new"]
            tbl_data = self.get_tbl_data(key)
            manual_data = self.get_manual_tbl_data(key)
            if tbl_data:
                for i, label_text in enumerate(labels):
                    l = tk.Label(self.transform_frame, text=f"{label_text}: {tbl_data[i]}", anchor="w")
                    l.grid(row=i, column=0, padx=5, pady=2, sticky="w")
            else:
                l = tk.Label(self.transform_frame, text="No transformation data found", anchor="w")
                l.grid(row=0, column=0, padx=5, pady=2, sticky="w")
            if manual_data:
                for i, label_text in enumerate(labels):
                    l = tk.Label(self.transform_frame, text=f"manual {label_text}: {manual_data[i]}", anchor="w")
                    l.grid(row=i, column=1, padx=15, pady=2, sticky="w")
            else:
                l = tk.Label(self.transform_frame, text="No manual transformation data", anchor="w")
                l.grid(row=0, column=1, padx=15, pady=2, sticky="w")
        else:
            self.master.title("No images found")
            self.file_info_label.config(text="")
            self.polarity_label.config(text="")
            for widget in self.transform_frame.winfo_children():
                widget.destroy()

        if self.object_viewer_window is not None and tk.Toplevel.winfo_exists(self.object_viewer_window):
            self.update_object_viewer()

    def next_image(self, event):
        if images:
            self.current_index = (self.current_index + 1) % len(images)
            self.show_image()

    def prev_image(self, event):
        if images:
            self.current_index = (self.current_index - 1) % len(images)
            self.show_image()

    def delete_manual(self):
        if not images:
            return
        image_path = images[self.current_index]
        file_name = os.path.basename(image_path)
        key = self.extract_key(file_name)
        manual_path = os.path.join("particles", key, "manual_xform.tbl")
        if os.path.exists(manual_path):
            try:
                os.remove(manual_path)
                print(f"Deleted {manual_path}")
            except Exception as e:
                print(f"Error deleting manual_xform.tbl: {e}")
        else:
            print("No manual file to delete.")
        self.show_image()

    def manual_correct(self):
        if not images:
            return
        image_path = images[self.current_index]
        file_name = os.path.basename(image_path)
        key = self.extract_key(file_name)
        particles_dir = os.path.join("particles", key)
        manual_path = os.path.join(particles_dir, "manual_xform.tbl")
        if os.path.exists(manual_path):
            try:
                with open(manual_path, "r") as f:
                    data_line = f.readline().strip()
            except Exception as e:
                print(f"Error reading manual_xform.tbl: {e}")
                return
        else:
            data_line = "0 0 0 0 0 0"
        parts = data_line.split()
        if len(parts) < 6:
            parts = ["0"] * 6

        editor = tk.Toplevel(self.master)
        editor.title(f"Manual Correction for {key}")
        menubar = tk.Menu(editor)
        view_menu = tk.Menu(menubar, tearoff=0)
        view_menu.add_command(label="Increase Image Size (+)", command=self.increase_image_size)
        view_menu.add_command(label="Decrease Image Size (-)", command=self.decrease_image_size)
        view_menu.add_command(label="Default Image Size (0)", command=self.default_image_size)
        view_menu.add_separator()
        menubar.add_cascade(label="View", menu=view_menu)
        editor.config(menu=menubar)
        # Bind keys in the manual correction window
        editor.bind("<KeyPress-equal>", lambda event: self.increase_image_size())
        editor.bind("<KeyPress-KP_Add>", lambda event: self.increase_image_size())
        editor.bind("<KeyPress-minus>", lambda event: self.decrease_image_size())
        editor.bind("<KeyPress-0>", lambda event: self.default_image_size())
        # The Reference Settings dialog will be opened via the View menu command below.
        # (update_canvas_image will be defined shortly)
        
        view_menu.add_command(label="Reference Settings...", command=lambda: self.open_ref_settings_dialog(editor, update_canvas_image))
        
        try:
            init_dx = float(parts[0])
            init_dy = float(parts[1])
            init_dtheta = float(parts[5])
            init_flip = int(parts[4])
        except Exception:
            init_dx, init_dy, init_dtheta, init_flip = 0.0, 0.0, 0.0, 0
        if init_flip == 180:
            init_dtheta = -init_dtheta
        transform = {"dx": init_dx, "dy": init_dy, "dtheta": init_dtheta}
        flip_value = init_flip
        try:
            original_image = Image.open(image_path)
        except Exception as e:
            print(f"Error opening image {image_path}: {e}")
            editor.destroy()
            return

        canvas = tk.Canvas(editor, width=350, height=350, bg="gray")
        canvas.grid(row=1, column=0, columnspan=3, padx=5, pady=5)
        translation_label = tk.Label(editor, text="")
        translation_label.grid(row=2, column=0, columnspan=3)
        rotation_label = tk.Label(editor, text="")
        rotation_label.grid(row=3, column=0, columnspan=3)
        display_image = {}

        def update_manual_entries():
            dx_entry.delete(0, tk.END)
            dx_entry.insert(0, f"{transform['dx']:.2f}")
            dy_entry.delete(0, tk.END)
            dy_entry.insert(0, f"{transform['dy']:.2f}")
            dtheta_entry.delete(0, tk.END)
            dtheta_entry.insert(0, f"{transform['dtheta']:.2f}")

        def update_canvas_image():
            abs_rot = transform["dtheta"]
            rotated = original_image.rotate(-abs_rot, expand=True).convert("RGBA")
            if self.scale_factor != 1.0:
                new_size = (int(rotated.width * self.scale_factor), int(rotated.height * self.scale_factor))
                rotated = rotated.resize(new_size, Image.ANTIALIAS)
            if flip_value == 180:
                rotated = rotated.transpose(Image.FLIP_TOP_BOTTOM)
            buffer = io.BytesIO()
            rotated.save(buffer, format="PNG")
            encoded = base64.b64encode(buffer.getvalue())
            tk_rotated = tk.PhotoImage(data=encoded.decode("ascii"))
            display_image["rotated"] = tk_rotated
            margin = 50
            img_w, img_h = rotated.size
            canvas_width = img_w + margin
            canvas_height = img_h + margin
            canvas.config(width=canvas_width, height=canvas_height)
            canvas.delete("all")
            center_x = canvas_width // 2
            center_y = canvas_height // 2
            canvas.create_image(
                center_x + transform["dx"] * self.scale_factor,
                center_y - transform["dy"] * self.scale_factor,
                image=tk_rotated
            )
            try:
                # Load and process the reference image with opacity and cut off thresholds
                ref_img = Image.open(REFERENCE_IMAGE_PATH).convert("RGBA")
                if self.scale_factor != 1.0:
                    new_size = (int(ref_img.width * self.scale_factor), int(ref_img.height * self.scale_factor))
                    ref_img = ref_img.resize(new_size, Image.ANTIALIAS)
                # Create a mask: only pixels with gray value between ref_black and ref_white get self.ref_alpha
                gray = ref_img.convert("L")
                def threshold_func(x):
                    return self.ref_alpha if (self.ref_black <= x <= self.ref_white) else 0
                mask = gray.point(threshold_func)
                r, g, b, _ = ref_img.split()
                ref_img = Image.merge("RGBA", (r, g, b, mask))
                buffer_ref = io.BytesIO()
                ref_img.save(buffer_ref, format="PNG")
                encoded_ref = base64.b64encode(buffer_ref.getvalue())
                tk_ref = tk.PhotoImage(data=encoded_ref.decode("ascii"))
                display_image["ref"] = tk_ref
                canvas.create_image(center_x, center_y, image=tk_ref)
            except Exception:
                pass
            translation_label.config(
                text=f"Translation delta: dx = {transform['dx']:.2f}, dy = {transform['dy']:.2f}"
            )
            rotation_label.config(text=f"Rotation delta: dθ = {transform['dtheta']:.2f}")
            update_manual_entries()

        drag_data = {"x": 0, "y": 0}
        def on_button_press(event):
            drag_data["x"] = event.x
            drag_data["y"] = event.y
        def on_mouse_drag(event):
            delta_x = event.x - drag_data["x"]
            delta_y = event.y - drag_data["y"]
            transform["dx"] += delta_x / self.scale_factor
            transform["dy"] -= delta_y / self.scale_factor
            drag_data["x"] = event.x
            drag_data["y"] = event.y
            update_canvas_image()
        canvas.bind("<ButtonPress-1>", on_button_press)
        canvas.bind("<B1-Motion>", on_mouse_drag)
        def update_rotation(val):
            try:
                transform["dtheta"] = float(val)
            except Exception:
                transform["dtheta"] = 0
            update_canvas_image()
        rotation_scale = tk.Scale(
            editor, from_=-180, to=180, orient="horizontal",
            label="Rotation Delta", command=update_rotation
        )
        rotation_scale.set(transform["dtheta"])
        rotation_scale.grid(row=4, column=0, columnspan=3, sticky="ew", padx=5, pady=5)
        manual_frame = tk.Frame(editor)
        manual_frame.grid(row=5, column=0, columnspan=3, sticky="ew", padx=5, pady=5)
        dy_text = "dy:" if flip_value != 180 else "dy:"
        dtheta_text = "dθ:" if flip_value != 180 else "-dθ:"
        tk.Label(manual_frame, text="dx:").grid(row=0, column=0, padx=2)
        dx_entry = tk.Entry(manual_frame, width=10)
        dx_entry.grid(row=0, column=1, padx=2)
        dy_label = tk.Label(manual_frame, text=dy_text)
        dy_label.grid(row=0, column=2, padx=2)
        dy_entry = tk.Entry(manual_frame, width=10)
        dy_entry.grid(row=0, column=3, padx=2)
        dtheta_label = tk.Label(manual_frame, text=dtheta_text)
        dtheta_label.grid(row=0, column=4, padx=2)
        dtheta_entry = tk.Entry(manual_frame, width=10)
        dtheta_entry.grid(row=0, column=5, padx=2)
        def on_dx_entry(event):
            try:
                transform["dx"] = float(dx_entry.get())
            except ValueError:
                pass
            update_canvas_image()
        def on_dy_entry(event):
            try:
                transform["dy"] = float(dy_entry.get())
            except ValueError:
                pass
            update_canvas_image()
        def on_dtheta_entry(event):
            try:
                transform["dtheta"] = float(dtheta_entry.get())
                rotation_scale.set(transform["dtheta"])
            except ValueError:
                pass
            update_canvas_image()
        dx_entry.bind("<Return>", on_dx_entry)
        dy_entry.bind("<Return>", on_dy_entry)
        dtheta_entry.bind("<Return>", on_dtheta_entry)
        flip_frame = tk.Frame(editor)
        flip_frame.grid(row=6, column=0, columnspan=3, pady=5)
        flip_label = tk.Label(flip_frame, text=f"x-rot: {flip_value}")
        flip_label.pack(side="left", padx=5)
        def flip_action():
            nonlocal flip_value
            flip_value = 180 if flip_value == 0 else 0
            flip_label.config(text=f"x-rot: {flip_value}")
            dtheta_label.config(text="-dθ:" if flip_value == 180 else "dθ:")
            update_canvas_image()
        flip_button = tk.Button(flip_frame, text="Flip", command=flip_action)
        flip_button.pack(side="left", padx=5)
        def save_manual():
            manual_path = os.path.join("particles", key, "manual_xform.tbl")
            if transform["dx"] == 0 and transform["dy"] == 0 and transform["dtheta"] == 0 and flip_value == 0:
                if os.path.exists(manual_path):
                    try:
                        os.remove(manual_path)
                    except Exception as e:
                        print(f"Error deleting manual_xform.tbl: {e}")
                editor.destroy()
                if editor in self.manual_editors:
                    self.manual_editors.remove(editor)
                self.show_image()
                return
            if flip_value == 180:
                dy_to_save = transform["dy"]
                dtheta_to_save = -transform["dtheta"]
            else:
                dy_to_save = transform["dy"]
                dtheta_to_save = transform["dtheta"]
            new_values = [
                str(transform["dx"]),
                str(dy_to_save),
                parts[2],
                parts[3],
                str(flip_value),
                str(dtheta_to_save)
            ]
            try:
                with open(manual_path, "w") as f:
                    f.write(" ".join(new_values) + "\n")
            except Exception as e:
                print(f"Error saving manual_xform.tbl: {e}")
            editor.destroy()
            if editor in self.manual_editors:
                self.manual_editors.remove(editor)
            self.show_image()
        save_button = tk.Button(editor, text="Save", command=save_manual)
        save_button.grid(row=7, column=0, columnspan=3, pady=10)
        def on_close():
            manual_path = os.path.join("particles", key, "manual_xform.tbl")
            if transform["dx"] == 0 and transform["dy"] == 0 and transform["dtheta"] == 0 and flip_value == 0:
                if os.path.exists(manual_path):
                    try:
                        os.remove(manual_path)
                    except Exception as e:
                        print(f"Error deleting manual_xform.tbl: {e}")
            editor.destroy()
            if editor in self.manual_editors:
                self.manual_editors.remove(editor)
            self.show_image()
        editor.protocol("WM_DELETE_WINDOW", on_close)
        update_canvas_image()
        editor.update_canvas = update_canvas_image
        self.manual_editors.append(editor)

    def open_ref_settings_dialog(self, parent, update_callback):
        dialog = tk.Toplevel(parent)
        dialog.title("Reference Settings")
        tk.Label(dialog, text="Reference Opacity (%)").pack(padx=10, pady=5)
        opacity_scale = tk.Scale(dialog, from_=0, to=100, orient="horizontal", resolution=1,
                                 command=lambda val: self.set_ref_opacity(val, update_callback))
        opacity_scale.set(int(self.ref_alpha / 255 * 100))
        opacity_scale.pack(padx=10, pady=5, fill="x")
        
        tk.Label(dialog, text="Black Cut Off").pack(padx=10, pady=5)
        black_scale = tk.Scale(dialog, from_=0, to=255, orient="horizontal", resolution=1,
                               command=lambda val: self.set_ref_black(val, update_callback))
        black_scale.set(self.ref_black)
        black_scale.pack(padx=10, pady=5, fill="x")
        
        tk.Label(dialog, text="White Cut Off").pack(padx=10, pady=5)
        white_scale = tk.Scale(dialog, from_=0, to=255, orient="horizontal", resolution=1,
                               command=lambda val: self.set_ref_white(val, update_callback))
        white_scale.set(self.ref_white)
        white_scale.pack(padx=10, pady=5, fill="x")

    def set_ref_opacity(self, val, update_callback):
        self.ref_alpha = int(float(val) / 100 * 255)
        update_callback()

    def set_ref_black(self, val, update_callback):
        self.ref_black = int(val)
        update_callback()

    def set_ref_white(self, val, update_callback):
        self.ref_white = int(val)
        update_callback()

    def increase_image_size(self):
        self.scale_factor *= 1.1
        self.show_image()
        self.update_manual_editors_scale()

    def decrease_image_size(self):
        self.scale_factor /= 1.1
        self.show_image()
        self.update_manual_editors_scale()

    def default_image_size(self):
        self.scale_factor = 1.0
        self.show_image()
        self.update_manual_editors_scale()

    def update_manual_editors_scale(self):
        for editor in self.manual_editors:
            try:
                editor.update_canvas()
            except Exception:
                pass

    def open_go_window(self):
        go_win = tk.Toplevel(self.master)
        go_win.title("Go to Image")
        listbox = tk.Listbox(go_win, width=50, height=15)
        scrollbar = tk.Scrollbar(go_win, orient="vertical", command=listbox.yview)
        listbox.config(yscrollcommand=scrollbar.set)
        listbox.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        for i, path in enumerate(images):
            listbox.insert("end", os.path.basename(path))
        def on_select(event):
            selection = listbox.curselection()
            if selection:
                index = selection[0]
                self.current_index = index
                self.show_image()
                go_win.destroy()
        listbox.bind("<Double-Button-1>", on_select)

    def go_to_image(self, path):
        try:
            idx = images.index(path)
            self.current_index = idx
            self.show_image()
        except ValueError:
            pass

    def object_viewer(self):
        current_path = images[self.current_index]
        current_file = os.path.basename(current_path)
        current_key = self.extract_key(current_file)
        family = current_key.rsplit('_', 1)[0]
        if self.object_viewer_window is None or not tk.Toplevel.winfo_exists(self.object_viewer_window):
            self.object_viewer_window = tk.Toplevel(self.master)
            self.object_viewer_window.title(f"Object Viewer: {family}")
            self.ov_current_family = family
            ov_menubar = tk.Menu(self.object_viewer_window)
            ov_view_menu = tk.Menu(ov_menubar, tearoff=0)
            ov_view_menu.add_command(label="Increase Image Size (+)", command=self.ov_increase_image_size)
            ov_view_menu.add_command(label="Decrease Image Size (-)", command=self.ov_decrease_image_size)
            ov_view_menu.add_command(label="Default Image Size (0)", command=self.ov_default_image_size)
            ov_menubar.add_cascade(label="View", menu=ov_view_menu)
            self.object_viewer_window.config(menu=ov_menubar)
            self.object_viewer_window.bind("<KeyPress-equal>", lambda event: self.ov_increase_image_size())
            self.object_viewer_window.bind("<KeyPress-KP_Add>", lambda event: self.ov_increase_image_size())
            self.object_viewer_window.bind("<KeyPress-minus>", lambda event: self.ov_decrease_image_size())
            self.object_viewer_window.bind("<KeyPress-0>", lambda event: self.ov_default_image_size())
            self.ov_frame = tk.Frame(self.object_viewer_window)
            self.ov_frame.pack(fill="both", expand=True)
            self.update_object_viewer(force=True)
        else:
            if self.ov_current_family != family:
                self.ov_current_family = family
                self.update_object_viewer(force=True)

    def update_object_viewer(self, force=False):
        if self.object_viewer_window is None or not tk.Toplevel.winfo_exists(self.object_viewer_window):
            self.object_viewer_window = None
            self.ov_current_family = None
            return
        current_path = images[self.current_index]
        current_file = os.path.basename(current_path)
        current_key = self.extract_key(current_file)
        family = current_key.rsplit('_', 1)[0]
        if not force and self.ov_current_family == family:
            return
        self.ov_current_family = family
        self.object_viewer_window.title(f"Object Viewer: {family}")
        for widget in self.ov_frame.winfo_children():
            widget.destroy()
        family_files = []
        for path in images:
            file_name = os.path.basename(path)
            key = self.extract_key(file_name)
            if key.startswith(family + "_"):
                family_files.append((key, path))
        family_files.sort()
        columns = 3
        for i, (key, path) in enumerate(family_files):
            obj_frame = tk.Frame(self.ov_frame, bd=2, relief="groove")
            obj_frame.grid(row=i // columns, column=i % columns, padx=5, pady=5)
            name_label = tk.Label(obj_frame, text=key)
            name_label.pack()
            polarity_val = polarity_map.get(key, "N/A")
            fg_color = "black"
            if polarity_val == "1":
                fg_color = COLOR_1
            elif polarity_val == "0":
                fg_color = COLOR_0
            polarity_label = tk.Label(obj_frame, text=f"Polarity: {polarity_val}", fg=fg_color)
            polarity_label.pack()
            img = self.load_and_scale_image_OV(path)
            image_label = tk.Label(obj_frame, image=img)
            image_label.image = img
            image_label.pack()
            obj_frame.bind("<Button-1>", lambda event, p=path: self.go_to_image(p))
            for child in obj_frame.winfo_children():
                child.bind("<Button-1>", lambda event, p=path: self.go_to_image(p))

    def ov_increase_image_size(self):
        self.ov_scale_factor *= 1.1
        self.update_object_viewer(force=True)

    def ov_decrease_image_size(self):
        self.ov_scale_factor /= 1.1
        self.update_object_viewer(force=True)

    def ov_default_image_size(self):
        self.ov_scale_factor = 1.0
        self.update_object_viewer(force=True)


if __name__ == '__main__':
    root = tk.Tk()
    app = PhotoBrowser(root)
    root.mainloop()

