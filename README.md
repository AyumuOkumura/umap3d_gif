# umap3dGIF 🎬

Create stunning rotating 3D UMAP visualizations as animated GIFs from Seurat single-cell RNA-seq objects.

[![License:  MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

`umap3dGIF` is an R package that generates high-quality, rotating 3D UMAP plots as animated GIFs from Seurat objects. It's designed to be memory-efficient, customizable, and easy to use for showcasing single-cell clustering results.

### Features

- ✨ **Automatic 3D UMAP calculation** or use existing reductions
- 🎨 **Flexible coloring options** with support for scCustomize palettes
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
