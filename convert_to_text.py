import tkinter as tk
from tkinter import filedialog
import os
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
from tqdm import tqdm
import concurrent.futures
from PIL import Image, ImageTk
from tkinter.ttk import Progressbar
from tkinter import Scale
from concurrent.futures import ThreadPoolExecutor

# Define a global variable to store the paths of extracted text files
extracted_files = []

# Existing code for GUI setup, function definitions, etc.

def choose_folder():
    root = tk.Tk()
    root.withdraw()  # Hide the main window
    folder_path = filedialog.askdirectory()  # Open a dialog to select a folder
    return folder_path

def plot_selected_files():
    # Create a loading message window
    loading_message = tk.Toplevel(root)
    loading_message.title("Loading...")
    loading_message.geometry("200x100")
    loading_label = tk.Label(loading_message, text="Loading, please wait...")
    loading_label.pack()

    # Force the window to update and display
    loading_message.update()

    # Get the selected files from the listbox
    selected_file_indices = listbox.curselection()
    if selected_file_indices:
        # Clear the plot before plotting the new files
        plt.clf()

        # Initialize arrays to store longitude, latitude, and SSHA values
        merged_lons = []
        merged_lats = []
        merged_ssha = []

        # Function to plot each selected file
        def plot_file(selected_file_info):
            nonlocal merged_lons, merged_lats, merged_ssha

            selected_file_path = selected_file_info[0]
            print("Plotting selected file:", selected_file_path)

            # Read the NetCDF file using xarray
            ds_expert = xr.open_dataset(selected_file_path)
            ds_expert = ds_expert.assign_coords(longitude=(((ds_expert.longitude + 180) % 360) - 180))

            # Extract longitude, latitude, and SSHA variables
            lon = ds_expert['longitude'].values
            lat = ds_expert['latitude'].values
            ssha = ds_expert['ssha'].values

            # Append longitude, latitude, and SSHA values to the merged arrays
            merged_lons.extend(lon)
            merged_lats.extend(lat)
            merged_ssha.extend(ssha)

        # Get selected file info for all selected files
        selected_files_info = [listbox.get(index) for index in selected_file_indices]

        # Plot selected files in parallel using ThreadPoolExecutor
        with ThreadPoolExecutor() as executor:
            executor.map(plot_file, selected_files_info)

        # Convert merged lists to numpy arrays
        merged_lons = np.array(merged_lons)
        merged_lats = np.array(merged_lats)
        merged_ssha = np.array(merged_ssha)

        # Get latitude and longitude ranges from entry fields
        min_lat = float(min_lat_entry.get())
        max_lat = float(max_lat_entry.get())
        min_lon = float(min_lon_entry.get())
        max_lon = float(max_lon_entry.get())

        vmin = vmin_scale.get()
        vmax = vmax_scale.get()

    # Plot the merged SSHA values on a single world map with specified latitude and longitude ranges
        merged_fig = plot_data(merged_lons, merged_lats, merged_ssha, min_lat, max_lat, min_lon, max_lon, vmin, vmax)
        update_plot_canvas(merged_fig)   # Update the plot on the canvas

        # Close the loading message window after plotting is completed
        loading_message.destroy()
    else:
        print("No file selected.")

def filter_files_by_coordinates(folder_path, min_lat, max_lat, min_lon, max_lon, pbar=None):
    # List to store filtered file paths
    filtered_files = []

    # Total number of files
    all_files = len(os.listdir(folder_path))

    # Iterate through files in the folder
    for idx, file_name in enumerate(os.listdir(folder_path), 1):
        file_path = os.path.join(folder_path, file_name)

        # Update progress bar
        if pbar is not None:
            pbar['value'] = idx / all_files * 100
            pbar.update()

        # Check if the file is a NetCDF file
        if file_name.endswith('.nc'):
            try:
                # Open the NetCDF file
                with xr.open_dataset(file_path) as dataset:
                    # Extract latitude and longitude coordinates from the dataset
                    latitudes = dataset['latitude'].values
                    longitudes = dataset['longitude'].values
                    ssha = dataset['ssha'].values  # Assuming 'ssha' is the variable name for Sea Surface Height Anomaly

                    # Find the indices corresponding to latitude and longitude ranges
                    lat_indices = np.where((latitudes >= min_lat) & (latitudes <= max_lat))[0]
                    lon_indices = np.where((longitudes >= min_lon) & (longitudes <= max_lon))[0]

                    # Check for common indices (intersection) for latitude and longitude
                    common_indices = np.intersect1d(lat_indices, lon_indices)

                    # Filter the latitude, longitude, and SSHA values corresponding to common indices
                    filtered_latitudes = latitudes[common_indices]
                    filtered_longitudes = longitudes[common_indices]
                    filtered_ssha = ssha[common_indices]

                    # Check if there are non-NaN SSHA values within the specified latitude and longitude ranges
                    if not np.isnan(filtered_ssha).all():
                        filtered_files.append((file_path, filtered_latitudes, filtered_longitudes))
            except (OSError, RuntimeError) as e:
                print(f"Error processing file {file_path}: {e}")
                continue

    return filtered_files

def button_3_clicked():
    # Get latitude and longitude values from the entry fields
    min_lat = float(min_lat_entry.get())
    max_lat = float(max_lat_entry.get())
    min_lon = float(min_lon_entry.get())
    max_lon = float(max_lon_entry.get())

    # Filter files based on the provided coordinates
    folder_path = choose_folder()  # Choose the folder path through file management system

    # Remove previous progress bar if it exists
    for widget in root.winfo_children():
        if isinstance(widget, Progressbar):
            widget.destroy()

    # Create a progress bar
    pbar = Progressbar(root, orient="horizontal", length=200, mode="determinate")
    pbar.pack(side="bottom")

    filtered_files = filter_files_by_coordinates(folder_path, min_lat, max_lat, min_lon, max_lon, pbar)
    valid_files = len(filtered_files)

    print("Filtered files:")
    for file_name in filtered_files:
        print(file_name)
    print("Number of valid files:", valid_files)
    print("Total files:", len(os.listdir(folder_path)))

    # Update the listbox with the filtered files
    listbox.delete(0, tk.END)  # Clear the listbox
    for file_name in filtered_files:
        listbox.insert(tk.END, file_name)

    # Update the label with the count of valid files
    valid_files_label.config(text=f"Valid Files: {valid_files}/{len(os.listdir(folder_path))}")

def plot_data(lons, lats, ssha_values, min_lat, max_lat, min_lon, max_lon, vmin, vmax):
    # Define the path to the coastline shapefile
    coastline_shapefile = "/home/abcd/Downloads/ne_10m_coastline/ne_10m_coastline.shp"

    # Check if the shapefile exists
    if os.path.exists(coastline_shapefile):
        # Read the coastline shapefile
        coastlines = gpd.read_file(coastline_shapefile)

        # Create a figure and axis using Cartopy
        fig, ax = plt.subplots(figsize=(12, 6), subplot_kw=dict(projection=ccrs.PlateCarree()))

        # Add blue marble background
        ax.stock_img()

        # Plot SSHA values within the specified latitude and longitude range
        pcm = ax.scatter(lons, lats, c=ssha_values, cmap='Spectral_r', vmin=vmin, vmax=vmax, s=5, transform=ccrs.PlateCarree())
        ax.set_extent([min_lon, max_lon, min_lat, max_lat], crs=ccrs.PlateCarree())

        # Add coastlines
        coastlines.plot(ax=ax, color='black', linewidth=0.5)

        # Add colorbar
        cbar = plt.colorbar(pcm, ax=ax, orientation='vertical', shrink=0.7, label='SSHA [units]')

        # Add gridlines
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False

        # Add title and labels
        ax.set_title('Sea Surface Height Anomalies')
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')

        # Return the figure
        return fig
    else:
        print("Coastline shapefile not found.")

def update_plot_canvas(fig):
    # Clear the previous plot on canvas
    for widget in center_frame.winfo_children():
        widget.destroy()

    # Create a new matplotlib canvas for embedding the updated plot
    canvas_new = FigureCanvasTkAgg(fig, master=center_frame)
    canvas_new.draw()
    canvas_new.get_tk_widget().pack(expand=True, fill=tk.BOTH)

    # Add Matplotlib navigation toolbar
    toolbar = NavigationToolbar2Tk(canvas_new, center_frame)
    toolbar.update()
    canvas_new._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

def choose_files():
    root = tk.Tk()
    root.withdraw()
    file_paths = filedialog.askopenfilenames(filetypes=[("NetCDF files", "*.nc")])
    root.destroy()
    return file_paths

def extract_data_to_text(nc_file_path, output_folder, min_lat, max_lat, min_lon, max_lon):
    try:
        # Open the NetCDF file using xarray
        ds = xr.open_dataset(nc_file_path)

        # Extract latitude, longitude, and SSHA data
        latitudes = ds['latitude'].values
        longitudes = ds['longitude'].values
        ssha_values = ds['ssha'].values

        # Flatten arrays and get corresponding indices of non-NaN values in SSHA
        valid_indices = np.where(~np.isnan(ssha_values))

        # Filter data based on user-defined latitude and longitude ranges
        latitudes = latitudes[valid_indices].flatten()  # Ensure 1-dimensional
        longitudes = longitudes[valid_indices].flatten()  # Ensure 1-dimensional
        ssha_values = ssha_values[valid_indices].flatten()  # Ensure 1-dimensional

        latitudes = np.round(latitudes, 2)  # Round latitude values to two decimal places
        longitudes = np.round(longitudes, 2)  # Round longitude values to two decimal places

        # Construct the output file path
        output_file_path = os.path.join(output_folder, os.path.basename(nc_file_path).replace('.nc', '_extracted_data.txt'))

        # Write latitude, longitude, and SSHA data to the text file
        with open(output_file_path, 'w') as file:
            file.write("Latitude\tLongitude\tSSHA\n")
            for lat, lon, ssha in zip(latitudes, longitudes, ssha_values):
                file.write(f"{lat}\t{lon}\t{ssha}\n")

        print(f"Extracted data from {nc_file_path} saved to: {output_file_path}")

        # Return the file path
        return output_file_path

    except Exception as e:
        print(f"Error extracting data from {nc_file_path}: {e}")
        return ''  # Return an empty string if an error occurs



def save_selected_files_as_text():
    # Get the indices of the selected files in the listbox
    selected_indices = listbox.curselection()
   
    if selected_indices:
        # Create a folder to save the text files if it doesn't exist
        output_folder = filedialog.askdirectory(title="Select output folder")
        if output_folder:
            # Prompt the user for the custom file name
            custom_file_name = tk.simpledialog.askstring("Custom File Name", "Enter a custom file name:")
            if custom_file_name:
                combined_file_name = f"{custom_file_name}.txt"
            else:
                combined_file_name = "combined_files.txt"
           
            # Initialize a list to store content of all text files
            all_files_content = []
           
            # Get latitude and longitude values from the entry fields
            min_lat = float(min_lat_entry.get())
            max_lat = float(max_lat_entry.get())
            min_lon = float(min_lon_entry.get())
            max_lon = float(max_lon_entry.get())
           
            for index in selected_indices:
                file_info = listbox.get(index)
                nc_file_path = file_info[0]
                file_content = extract_data_to_text(nc_file_path, output_folder, min_lat, max_lat, min_lon, max_lon)
                all_files_content.append(file_content)
           
            # Write the content of all text files into a single file
            combined_file_path = os.path.join(output_folder, combined_file_name)
            with open(combined_file_path, "w") as combined_file:
                for file_content in all_files_content:
                    with open(file_content, 'r') as file:
                        combined_file.write(file.read())
           
            print(f"Text files saved and combined successfully as: {combined_file_path}")
    else:
        print("No file selected.")

        
# Create the main window
root = tk.Tk()
root.title("SWOT_GUI")
root.minsize(1660, 1000)

# Load and resize the first image (adjust path and size as needed)
image_path1 = "/home/abcd/Downloads/download.jpeg"
image1 = Image.open(image_path1).resize((60, 60))  # Adjust the resize dimensions
photo1 = ImageTk.PhotoImage(image1)

# Load and resize the second image (adjust path and size as needed)
image_path2 = "/home/abcd/Downloads/Indian_Space_Research_Organisation_Logo.svg.png"
image2 = Image.open(image_path2).resize((60, 60))  # Adjust the resize dimensions
photo2 = ImageTk.PhotoImage(image2)

# Create top frame with sky blue background
top_frame = tk.Frame(root, width=200, height=40, pady=10, bg="sky blue")
top_frame.pack(side="top", fill="both")

# Add the first image to the top-left corner of the top frame
image_label1 = tk.Label(top_frame, image=photo1)
image_label1.image = photo1  # Store a reference to the image
image_label1.pack(side="left", anchor="nw")

# Add the second image next to the first image in the top frame
image_label2 = tk.Label(top_frame, image=photo2)
image_label2.image = photo2  # Store a reference to the image
image_label2.pack(side="left", anchor="nw")

# # Create left frame with sky blue background
# left_frame = tk.Frame(root, width=500, height=40, padx=30, bg="sky blue")
# left_frame.pack(side="right", fill="both")

# Load and display the GIF image in the left frame
gif_path = "/home/abcd/Downloads/karin_animation_copie.gif"  # Update the path to your GIF image
gif_image = Image.open(gif_path)

# Resize the GIF image to a smaller dimension
# Resize the GIF image to a smaller dimension
new_width = 50  # Adjust the width as needed
new_height = 50  # Adjust the height as needed
gif_resized = gif_image.resize((new_width, new_height), Image.BICUBIC)  # Use Image.BICUBIC as an alternative

# Convert the resized image to a Tkinter PhotoImage
gif_photo = ImageTk.PhotoImage(gif_resized)

# Calculate the width of the left frame
left_frame_width = new_width + 60  # Add extra padding as needed

# Create the left frame with the adjusted width
left_frame = tk.Frame(root, width=left_frame_width, height=new_height + 60, padx=30, bg="sky blue")
left_frame.pack(side="right", fill="both")

# Create a label to display the resized GIF image
gif_label = tk.Label(left_frame, image=gif_photo, bg="sky blue")
gif_label.pack(side="top", pady=20)

def update_image(ind):
    frame = gif_image.seek(ind)
    gif_frame = ImageTk.PhotoImage(gif_image)
    gif_label.configure(image=gif_frame)
    gif_label.image = gif_frame
    ind += 1
    root.after(100, update_image, ind % gif_image.n_frames)

# Start the GIF animation
update_image(0)

# Create fields for latitude and longitude with bold font
font_bold = ("Arial", 10, "bold")

# First line: Latitude labels and entries
latitude_frame = tk.Frame(left_frame, bg="sky blue")
latitude_frame.pack(side="top", fill="x", padx=5, pady=5)

min_lat_label = tk.Label(latitude_frame, text="Min Latitude:", font=font_bold)
min_lat_entry = tk.Entry(latitude_frame, font=font_bold)
max_lat_label = tk.Label(latitude_frame, text="Max Latitude:", font=font_bold)
max_lat_entry = tk.Entry(latitude_frame, font=font_bold)

min_lat_label.pack(side="left", padx=40, pady=5)
min_lat_entry.pack(side="left", padx=5, pady=5)
max_lat_label.pack(side="left", padx=10, pady=5)
max_lat_entry.pack(side="left", padx=5, pady=5)

# Second line: Longitude labels and entries
longitude_frame = tk.Frame(left_frame, bg="sky blue")
longitude_frame.pack(side="top", fill="x", padx=5, pady=5)

min_lon_label = tk.Label(longitude_frame, text="Min Longitude:", font=font_bold)
min_lon_entry = tk.Entry(longitude_frame, font=font_bold)
max_lon_label = tk.Label(longitude_frame, text="Max Longitude:", font=font_bold)
max_lon_entry = tk.Entry(longitude_frame, font=font_bold)

min_lon_label.pack(side="left", padx=40, pady=5)
min_lon_entry.pack(side="left", padx=5, pady=5)
max_lon_label.pack(side="left", padx=10, pady=5)
max_lon_entry.pack(side="left", padx=5, pady=5)

# Centering latitude and longitude frames
latitude_frame.pack_configure(anchor="center")
longitude_frame.pack_configure(anchor="center")

# Create buttons in left frame with bold font and sky blue color
button_font = ("Arial", 10, "bold")
button_bg_color = "sky blue"


# Create a frame to contain the buttons
button_frame = tk.Frame(left_frame, bg="sky blue")
button_frame.pack(side="top", padx=15, pady=25)

button3 = tk.Button(button_frame, text="Button 3", command=button_3_clicked, font=button_font, bg=button_bg_color)
button3.grid(row=0, column=0, padx=(10, 30), pady=5)  # Increased padding on the right

plot_ssha_button = tk.Button(button_frame, text="Plot SSHA", command=plot_selected_files, font=button_font, bg=button_bg_color)
plot_ssha_button.grid(row=0, column=1, padx=(30, 10), pady=5)  # Increased padding on the left

# Create a frame to contain the sliders
slider_frame = tk.Frame(left_frame, bg="sky blue")
slider_frame.pack(side="top", padx=15, pady=25)

# Add scales for adjusting vmin and vmax values
vmin_scale = Scale(slider_frame, from_=0, to=-0.5, resolution=0.01, orient=tk.HORIZONTAL, label="Min SSHA", font=font_bold)
vmin_scale.pack(side="left", padx=10, pady=5)

vmax_scale = Scale(slider_frame, from_=0, to=0.5, resolution=0.01, orient=tk.HORIZONTAL, label="Max SSHA", font=font_bold)
vmax_scale.pack(side="left", padx=10, pady=5)

# Create a frame to contain the "Save Selected as Text" button
save_button_frame = tk.Frame(left_frame, bg="sky blue")
save_button_frame.pack(side="top", padx=15, pady=25)

save_text_button = tk.Button(save_button_frame, text="Save Selected as Text", command=save_selected_files_as_text, font=button_font, bg=button_bg_color)
save_text_button.grid(row=0, column=0, padx=(30, 10), pady=5)

# Create a frame to contain the listbox
list_frame = tk.Frame(root, bg="white", width=200, height=200)
list_frame.pack(side="bottom", fill="both")

# Create a scrollbar for the listbox
scrollbar = tk.Scrollbar(list_frame, orient="vertical")

# Create a listbox widget to display the filtered files
listbox = tk.Listbox(list_frame, yscrollcommand=scrollbar.set, width=50, height=10, selectmode=tk.MULTIPLE)
scrollbar.config(command=listbox.yview)
scrollbar.pack(side="right", fill="y")
listbox.pack(side="left", fill="both", expand=True)

# Create a label to display the count of valid files
valid_files_label = tk.Label(root, text="Valid Files: 0/0")
valid_files_label.pack(side="bottom")

# Create center frame for displaying plots
center_frame = tk.Frame(root)
center_frame.pack(expand=True, fill="both")

# Load and display the image to be placed at the bottom center of the left frame
image_path_bottom = "/home/abcd/Downloads/output-onlinetexttools.png"  # Update the path to your image
image_bottom = Image.open(image_path_bottom)
image_bottom_resized = image_bottom.resize((200, 20))  # Adjust the size as needed
photo_bottom = ImageTk.PhotoImage(image_bottom_resized)

# Create a label for the image
image_label_bottom = tk.Label(left_frame, image=photo_bottom, bg="sky blue")
image_label_bottom.image = photo_bottom  # Store a reference to the image
image_label_bottom.pack(side="bottom", pady=20)  # Pack the label at the bottom center of the left frame


# Start the GUI event loop
root.mainloop()
