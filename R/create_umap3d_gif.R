#' Create 3D UMAP Rotating GIF (Memory Optimized Version)
#'
#' Generate a rotating 3D UMAP visualization as an animated GIF from a Seurat object.
#' This function creates a 3D UMAP projection (or uses an existing one) and captures
#' frames while rotating the view, then assembles them into a GIF with optional
#' background transparency.
#'
#' @param seurat_obj A Seurat object containing single-cell data
#' @param file_name Character.  Base name for the output GIF file (without extension).
#'   If NULL, uses the object name.  Default:  NULL
#' @param tag Character. Optional tag to append to the filename. Default: NULL
#' @param reduction Character. Name of the reduction to use as basis.  Default: "umap"
#' @param reduction. name. 3d Character. Name for the 3D reduction if calculating.  
#'   Default: "umap3d"
#' @param recalc_umap Logical. Force recalculation of 3D UMAP even if it exists.
#'   Default: FALSE
#' @param cluster_column Character. Name of metadata column containing cluster assignments.
#'   If NULL, uses Idents(seurat_obj). Default: NULL
#' @param palette Character. Color palette name from scCustomize.  Options include
#'   "varibow", "polychrome", "stepped", etc. Default: "varibow"
#' @param shuffle_colors Logical. Whether to shuffle the palette colors. Default: TRUE
#' @param dims Integer vector. Dimensions to use for UMAP calculation. Default: 1:30
#' @param n_components Integer. Number of UMAP components (must be 3 for 3D). Default: 3L
#' @param cols Character vector. Custom colors for clusters. If provided, overrides
#'   palette. Can be named (cluster names) or unnamed (by order). Default: NULL
#' @param label Logical. Whether to display cluster labels at centroids. Default: FALSE
#' @param label_size Numeric. Size of cluster labels. Default: 1.0
#' @param label_color Character. Color of cluster labels. Default: "gray"
#' @param fps Integer. Frames per second for the GIF.  Default: 10
#' @param duration Integer. Duration of the GIF in seconds. Default: 15
#' @param start_theta Numeric. Starting rotation angle (degrees). Default: 30
#' @param phi Numeric. Vertical viewing angle (degrees). Default: 30
#' @param zoom Numeric. Zoom level for the 3D view. Default: 0.6
#' @param point_size Numeric. Size of points in the 3D plot. Default: 3
#' @param window_size Integer vector of length 2. Width and height of the output.  
#'   Default: c(800, 800)
#' @param output_dir Character. Directory to save the GIF. Default: "."
#' @param bg_color Character. Background color for the plot. Default: "white"
#' @param transparency_fuzz Integer. Fuzz factor for background transparency (0-100).
#'   Higher values make more similar colors transparent. Set to 0 to disable.  Default: 10
#' @param cleanup Logical. Whether to delete temporary frame files after GIF creation.
#'   Default: TRUE
#' @param verbose Logical. Whether to print progress messages. Default: TRUE
#'
#' @return Invisibly returns the path to the generated GIF file.  
#'
#' @examples
#' \dontrun{
#' # Basic usage with default settings
#' create_umap3d_gif(seurat_obj)
#'
#' # Custom filename and cluster column
#' create_umap3d_gif(
#'   seurat_obj,
#'   file_name = "my_analysis",
#'   cluster_column = "cell_type"
#' )
#'
#' # With cluster labels and custom colors
#' create_umap3d_gif(
#'   seurat_obj,
#'   label = TRUE,
#'   cols = c("red", "blue", "green"),
#'   duration = 20,
#'   fps = 15
#' )
#'
#' # High quality with transparent background
#' create_umap3d_gif(
#'   seurat_obj,
#'   window_size = c(1200, 1200),
#'   bg_color = "white",
#'   transparency_fuzz = 15,
#'   point_size = 5
#' )
#' }
#'
#' @export
#' @importFrom Seurat RunUMAP Embeddings Idents
#' @importFrom rgl open3d par3d bg3d plot3d abclines3d text3d view3d rgl.snapshot close3d rgl.cur
#' @importFrom magick image_read image_transparent image_write
#' @importFrom gifski gifski
#' @importFrom scales hue_pal
create_umap3d_gif <- function(
    seurat_obj,
    file_name = NULL,
    tag = NULL,
    reduction = "umap",
    reduction.name.3d = "umap3d",
    recalc_umap = FALSE,
    cluster_column = NULL,
    palette = "varibow",
    shuffle_colors = TRUE,
    dims = 1:30,
    n_components = 3L,
    cols = NULL,
    label = FALSE,
    label_size = 1.0,
    label_color = "gray",
    fps = 10,
    duration = 15,
    start_theta = 30,
    phi = 30,
    zoom = 0.6,
    point_size = 3,
    window_size = c(800, 800),
    output_dir = ".",
    bg_color = "white",
    transparency_fuzz = 10,
    cleanup = TRUE,
    verbose = TRUE
) {
  
  # --- Input Validation ---
  validate_inputs(seurat_obj, fps, duration, n_components, window_size, 
                  transparency_fuzz, output_dir)
  
  # Setup cleanup handler
  on.exit({
    if (exists("temp_dir") && cleanup && dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
    }
    if (rgl:: rgl.cur() != 0) rgl::close3d()
  }, add = TRUE)
  
  # --- 1. Data Preparation ---
  target_reduction <- get_or_create_3d_reduction(
    seurat_obj, reduction, reduction.name.3d, recalc_umap, 
    dims, n_components, verbose
  )
  
  umap3d <- Seurat:: Embeddings(seurat_obj, target_reduction)[, 1:3]
  
  # --- 2. Cluster & Color Setup ---
  clusters <- get_cluster_assignments(seurat_obj, cluster_column)
  cell_colors <- assign_cell_colors(clusters, cols, palette, shuffle_colors, verbose)
  
  # --- 3. Plotting ---
  if (verbose) cat("Rendering 3D plot...\n")
  setup_3d_plot(umap3d, cell_colors, point_size, window_size, bg_color,
                clusters, label, label_color, label_size)
  
  # --- 4. Frame Capture (Memory Optimized) ---
  temp_dir <- create_temp_directory()
  capture_frames(temp_dir, fps, duration, start_theta, phi, zoom, 
                 bg_color, transparency_fuzz, verbose)
  
  # --- 5. GIF Assembly ---
  output_path <- assemble_gif(temp_dir, file_name, tag, output_dir, 
                               window_size, fps, verbose, seurat_obj)
  
  return(invisible(output_path))
}

# --- Helper Functions ---

#' Validate function inputs
#' @keywords internal
validate_inputs <- function(seurat_obj, fps, duration, n_components, 
                            window_size, transparency_fuzz, output_dir) {
  if (! inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required.  Install with: install.packages('Seurat')")
  }
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("rgl package is required. Install with: install.packages('rgl')")
  }
  if (!requireNamespace("magick", quietly = TRUE)) {
    stop("magick package is required. Install with: install.packages('magick')")
  }
  if (!requireNamespace("gifski", quietly = TRUE)) {
    stop("gifski package is required. Install with: install.packages('gifski')")
  }
  
  if (fps <= 0 || duration <= 0) {
    stop("fps and duration must be positive numbers")
  }
  if (n_components != 3L) {
    stop("n_components must be 3 for 3D visualization")
  }
  if (length(window_size) != 2 || any(window_size <= 0)) {
    stop("window_size must be a vector of two positive integers")
  }
  if (transparency_fuzz < 0 || transparency_fuzz > 100) {
    stop("transparency_fuzz must be between 0 and 100")
  }
  if (! dir.exists(output_dir)) {
    stop("output_dir does not exist:  ", output_dir)
  }
}

#' Get or create 3D reduction
#' @keywords internal
get_or_create_3d_reduction <- function(seurat_obj, reduction, reduction.name.3d,
                                        recalc_umap, dims, n_components, verbose) {
  if (! recalc_umap && reduction.name.3d %in% names(seurat_obj@reductions) && 
      ncol(seurat_obj[[reduction.name.3d]]) >= 3) {
    if (verbose) cat(sprintf("Using existing 3D reduction: %s\n", reduction.name.3d))
    return(reduction.name.3d)
  } else if (! recalc_umap && reduction %in% names(seurat_obj@reductions) && 
             ncol(seurat_obj[[reduction]]) >= 3) {
    if (verbose) cat(sprintf("Using existing 3D reduction: %s\n", reduction))
    return(reduction)
  } else {
    if (verbose) cat(sprintf("Calculating 3D UMAP into '%s'...\n", reduction.name.3d))
    seurat_obj <- Seurat::RunUMAP(
      seurat_obj, 
      dims = dims, 
      n.components = n_components,
      reduction.name = reduction.name.3d, 
      reduction.key = "UMAP3D_", 
      verbose = FALSE
    )
    return(reduction.name.3d)
  }
}

#' Get cluster assignments from Seurat object
#' @keywords internal
get_cluster_assignments <- function(seurat_obj, cluster_column) {
  if (is.null(cluster_column)) {
    return(factor(Seurat::Idents(seurat_obj)))
  } else {
    if (! cluster_column %in% colnames(seurat_obj@meta.data)) {
      stop("cluster_column '", cluster_column, "' not found in metadata")
    }
    return(factor(seurat_obj@meta.data[[cluster_column]]))
  }
}

#' Assign colors to cells
#' @keywords internal
assign_cell_colors <- function(clusters, cols, palette, shuffle_colors, verbose) {
  n_clusters <- length(levels(clusters))
  
  if (! is.null(cols)) {
    if (verbose) cat("Using provided custom colors\n")
    if (length(cols) < n_clusters) {
      warning("Number of colors provided is less than number of clusters.  Recycling colors.")
    }
    cell_colors <- if (! is.null(names(cols))) {
      cols[as.character(clusters)]
    } else {
      cols[as.numeric(clusters)]
    }
  } else if (requireNamespace("scCustomize", quietly = TRUE) && ! is.null(palette)) {
    if (verbose) cat(sprintf("Using scCustomize palette:  %s\n", palette))
    pal_colors <- scCustomize:: DiscretePalette_scCustomize(
      num_colors = n_clusters, 
      palette = palette, 
      shuffle_pal = shuffle_colors
    )
    cell_colors <- pal_colors[as.numeric(clusters)]
  } else {
    if (verbose && ! is.null(palette)) {
      cat("scCustomize not available, using default colors\n")
    }
    cell_colors <- scales::hue_pal()(n_clusters)[as.numeric(clusters)]
  }
  
  return(cell_colors)
}

#' Setup 3D plot
#' @keywords internal
setup_3d_plot <- function(umap3d, cell_colors, point_size, window_size, bg_color,
                          clusters, label, label_color, label_size) {
  options(rgl.useNULL = FALSE)
  rgl::open3d()
  rgl::par3d(windowRect = c(50, 50, 50 + window_size[1], 50 + window_size[2]))
  rgl::bg3d(bg_color)
  
  rgl::plot3d(umap3d, col = cell_colors, size = point_size, type = "p",
              axes = FALSE, box = FALSE, xlab = "", ylab = "", zlab = "")
  
  # Add coordinate axes
  rgl::abclines3d(0, 0, 0, a = 1, b = 0, c = 0, col = "gray", lwd = 1)
  rgl::abclines3d(0, 0, 0, a = 0, b = 1, c = 0, col = "gray", lwd = 1)
  rgl::abclines3d(0, 0, 0, a = 0, b = 0, c = 1, col = "gray", lwd = 1)
  
  # Add cluster labels if requested
  if (label) {
    add_cluster_labels(umap3d, clusters, label_color, label_size)
  }
}

#' Add cluster labels to 3D plot
#' @keywords internal
add_cluster_labels <- function(umap3d, clusters, label_color, label_size) {
  cluster_levels <- levels(clusters)
  for (lvl in cluster_levels) {
    idx <- which(clusters == lvl)
    if (length(idx) > 0) {
      centroid <- colMeans(umap3d[idx, , drop = FALSE])
      rgl::text3d(centroid[1], centroid[2], centroid[3], 
                  texts = lvl, col = label_color, cex = label_size, font = 2)
    }
  }
}

#' Create temporary directory for frames
#' @keywords internal
create_temp_directory <- function() {
  temp_dir <- file.path(tempdir(), paste0("umap_frames_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  return(temp_dir)
}

#' Capture frames for animation
#' @keywords internal
capture_frames <- function(temp_dir, fps, duration, start_theta, phi, zoom,
                           bg_color, transparency_fuzz, verbose) {
  n_frames <- fps * duration
  if (verbose) cat(sprintf("Generating %d frames in %s...\n", n_frames, temp_dir))
  
  thetas <- start_theta + 360 * (seq(1, n_frames) / n_frames)
  
  for (i in seq_len(n_frames)) {
    if (verbose && i %% 10 == 0) {
      cat(sprintf("\r  Progress: %d/%d", i, n_frames))
      flush.console()
    }
    
    rgl::view3d(theta = thetas[i], phi = phi, zoom = zoom)
    Sys.sleep(0.05) # GPU sync
    
    snap_file <- file.path(temp_dir, sprintf("frame_%04d.png", i))
    rgl::rgl.snapshot(snap_file)
    
    # Process transparency frame by frame (memory efficient)
    if (transparency_fuzz > 0) {
      process_frame_transparency(snap_file, bg_color, transparency_fuzz)
    }
    
    if (i %% 20 == 0) gc(verbose = FALSE)
  }
  
  if (verbose) cat("\nFrames captured.  Assembling GIF...\n")
}

#' Process frame transparency
#' @keywords internal
process_frame_transparency <- function(snap_file, bg_color, transparency_fuzz) {
  tryCatch({
    img <- magick::image_read(snap_file)
    img <- magick::image_transparent(img, color = bg_color, fuzz = transparency_fuzz)
    magick::image_write(img, snap_file)
    rm(img)
  }, error = function(e) {
    warning("Frame transparency processing failed: ", e$message)
  })
}

#' Assemble frames into GIF
#' @keywords internal
assemble_gif <- function(temp_dir, file_name, tag, output_dir, window_size, 
                         fps, verbose, seurat_obj) {
  obj_name <- sanitize_filename(file_name, seurat_obj, verbose)
  
  if (! is.null(tag) && nzchar(tag)) {
    obj_name <- paste(obj_name, sanitize_name(tag), sep = "_")
  }
  
  output_path <- file.path(output_dir, paste0(obj_name, "_umap3d.gif"))
  png_files <- list.files(temp_dir, full.names = TRUE, pattern = "\\.png$")
  
  if (length(png_files) == 0) {
    stop("No PNG frames were generated")
  }
  
  tryCatch({
    gifski:: gifski(
      png_files = sort(png_files),
      gif_file = output_path,
      width = window_size[1],
      height = window_size[2],
      delay = 1 / fps,
      loop = TRUE,
      progress = verbose
    )
    if (verbose) cat(sprintf("\nDone!  Saved to: %s\n", output_path))
  }, error = function(e) {
    stop("GIF generation failed: ", e$message)
  })
  
  return(output_path)
}

#' Sanitize filename
#' @keywords internal
sanitize_filename <- function(file_name, seurat_obj, verbose) {
  if (! is.null(file_name) && nzchar(file_name)) {
    return(sanitize_name(file_name))
  }
  
  obj_name_raw <- tryCatch(
    paste(deparse(substitute(seurat_obj)), collapse = ""),
    error = function(e) "seurat_object"
  )
  obj_name <- sanitize_name(obj_name_raw)
  
  if (nchar(obj_name) > 50 || obj_name == "") {
    if (verbose) message("Generated filename was invalid.  Using '3d_umap' as default.")
    obj_name <- "3d_umap"
  }
  
  return(obj_name)
}

#' Sanitize string for filename
#' @keywords internal
sanitize_name <- function(x) {
  gsub("[^[: alnum:]_.-]", "_", x)
}
