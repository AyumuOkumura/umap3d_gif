# umap3dGIF 🎬

Create stunning rotating 3D UMAP visualizations as animated GIFs from Seurat single-cell RNA-seq objects.

[![License:  MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

`umap3dGIF` is an R package that generates high-quality, rotating 3D UMAP plots as animated GIFs from Seurat objects. It's designed to be memory-efficient, customizable, and easy to use for showcasing single-cell clustering results.

### Features

- ✨ **Automatic 3D UMAP calculation** or use existing reductions
- 🎨 **Flexible coloring options** with customizable color schemes
- 🏷️ **Cluster labeling** with centroid positioning
- 🔄 **Smooth rotation animation** with customizable speed and duration
- 💾 **Memory-optimized** frame-by-frame processing
- 🎯 **Transparent backgrounds** for presentations
- 📦 **Export-ready** GIFs for publications and presentations

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("AyumuOkumura/umap3d_gif")
```

### System Dependencies

This package requires several dependencies:

```r
install.packages(c("Seurat", "rgl", "magick", "gifski", "scales"))
```

For `magick`, you may need to install ImageMagick on your system: 
- **Windows**: Download from [ImageMagick website](https://imagemagick.org/script/download.php)
- **macOS**: `brew install imagemagick`
- **Linux**: `sudo apt-get install libmagick++-dev`

## Quick Start

```r
library(umap3dGIF)
library(Seurat)

# Load your Seurat object
# seurat_obj <- readRDS("your_seurat_object.rds")

# Basic usage - creates a rotating 3D UMAP GIF
create_umap3d_gif(seurat_obj)
```

This will generate a file named `seurat_obj_umap3d.gif` in your current directory.

## Usage Examples

### Basic Usage with Custom Filename

```r
create_umap3d_gif(
  seurat_obj,
  file_name = "pbmc_analysis",
  output_dir = "./figures"
)
```

### With Cluster Labels

```r
create_umap3d_gif(
  seurat_obj,
  cluster_column = "cell_type",  # Use a specific metadata column
  label = TRUE,                  # Show cluster labels
  label_size = 1. 5,
  label_color = "black"
)
```

### Custom Colors

```r
# Define custom colors for your clusters
my_colors <- c("red", "blue", "green", "purple", "orange")

create_umap3d_gif(
  seurat_obj,
  cols = my_colors,
  point_size = 4
)
```

### High-Quality Output with Transparent Background

```r
create_umap3d_gif(
  seurat_obj,
  window_size = c(1200, 1200),  # Higher resolution
  bg_color = "white",
  transparency_fuzz = 15,       # Adjust for better transparency
  fps = 15,                     # Smoother animation
  duration = 20,                # Longer duration
  point_size = 5
)
```

### Advanced Customization

```r
create_umap3d_gif(
  seurat_obj,
  file_name = "my_cells",
  tag = "labeled",              # Adds tag to filename
  cluster_column = "seurat_clusters",
  recalc_umap = TRUE,          # Force recalculation of 3D UMAP
  dims = 1: 50,                 # Use more PCs
  label = TRUE,
  fps = 12,
  duration = 18,
  start_theta = 45,            # Starting rotation angle
  phi = 20,                    # Vertical viewing angle
  zoom = 0.7,
  window_size = c(1000, 1000),
  output_dir = "./output",
  cleanup = TRUE               # Remove temporary frames after creation
)
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `seurat_obj` | - | **Required. ** A Seurat object |
| `file_name` | `NULL` | Base name for output GIF (without extension) |
| `tag` | `NULL` | Optional tag to append to filename |
| `cluster_column` | `NULL` | Metadata column for clusters (uses `Idents()` if NULL) |
| `label` | `FALSE` | Show cluster labels at centroids |
| `cols` | `NULL` | Custom color vector for clusters |
| `fps` | `10` | Frames per second |
| `duration` | `15` | Animation duration in seconds |
| `window_size` | `c(800, 800)` | Output resolution (width, height) |
| `point_size` | `3` | Size of data points |
| `bg_color` | `"white"` | Background color |
| `transparency_fuzz` | `10` | Background transparency (0-100) |
| `recalc_umap` | `FALSE` | Force 3D UMAP recalculation |
| `dims` | `1:30` | Dimensions for UMAP calculation |
| `output_dir` | `"."` | Output directory |
| `verbose` | `TRUE` | Print progress messages |

## Output

The function creates a GIF file with the naming pattern:  `{file_name}_{tag}_umap3d.gif`

The GIF will show your 3D UMAP rotating 360 degrees, making it perfect for: 
- 📊 Presentations and talks
- 📄 Supplementary materials for publications
- 🌐 Lab websites and social media
- 📧 Sharing results with collaborators

## Tips for Best Results

1. **Resolution**: Use `window_size = c(1200, 1200)` or higher for publication-quality figures
2. **Smoothness**: Increase `fps` to 15-20 for smoother animations (larger file size)
3. **Duration**: 15-20 seconds is usually sufficient for one full rotation
4. **Transparency**:  Adjust `transparency_fuzz` if background removal isn't perfect (try values 5-20)
5. **Memory**: The function processes frames individually to minimize memory usage, but very high resolutions may still require substantial RAM

## Troubleshooting

### "rgl" window doesn't open
- Make sure you have a working display (X11 on Linux, XQuartz on macOS)
- On servers, you may need to set `options(rgl.useNULL = TRUE)`

### GIF looks choppy
- Increase `fps` (e.g., 15 or 20)
- Check that you have enough disk space for temporary frames

### Background isn't transparent
- Adjust `transparency_fuzz` (try values between 5-20)
- Ensure `bg_color` matches your actual plot background

### Out of memory errors
- Reduce `window_size`
- Decrease `duration` or `fps`
- Set `cleanup = TRUE` to remove temporary files immediately

## License

MIT License - see [LICENSE](LICENSE) file for details

## Citation

If you use this package in your research, please cite:

```
Okumura, A.  (2026). umap3dGIF: Create 3D UMAP Rotating GIF Animations from Seurat Objects. 
R package version 0.1.0. https://github.com/AyumuOkumura/umap3d_gif
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Issues

Found a bug or have a feature request?  Please open an issue on [GitHub](https://github.com/AyumuOkumura/umap3d_gif/issues).

## Acknowledgments

- Built on top of [Seurat](https://satijalab.org/seurat/)
- Uses [rgl](https://dmurdoch.github.io/rgl/) for 3D visualization
- GIF creation powered by [gifski](https://gif.ski/)
