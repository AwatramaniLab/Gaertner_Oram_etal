LoadVizgen_opt <- function(data.dir, 
                           fov = 'vz', 
                           assay = 'Vizgen',
                           mol.type = 'microns',
                           filter = '^Blank-',
                           z = 3L,
                           add.zIndex = TRUE, 
                           update.object = TRUE,
                           verbose,
                           ...)
{
  # reading data..
  data <- ReadVizgen_opt(data.dir = data.dir,
                         mol.type = mol.type,
                         filter = filter,
                         z = z,
                         verbose = verbose,
                         ...)
  
  if (verbose) { message("Creating Seurat object..") }  
  obj <- CreateSeuratObject(counts = data[["transcripts"]], assay = assay)
  
  # in case no segmentation is present, use boxes
  if (!"segmentations" %in% names(data)) {
    if ("boxes" %in% names(data)) {
      bound.boxes <- CreateSegmentation(data[["boxes"]])
      cents <- CreateCentroids(data[["centroids"]])
      bound.boxes.data <- list(centroids = cents, 
                               boxes = bound.boxes)
      if (verbose) { 
        message("Creating FOVs..", "\n",
                ">>> using box coordinates instead of segmentations") 
      }
      coords <- CreateFOV(coords = bound.boxes.data, 
                          type = c("boxes", "centroids"),
                          molecules = data[[mol.type]], 
                          assay = assay)
    } else { 
      # in case no segmentation & no boxes are present, use centroids only
      cents <- CreateCentroids(data[["centroids"]])
      if (verbose) { 
        message("Creating FOVs..", "\n",
                ">>> using only centroids") 
      }
      coords <- CreateFOV(coords = list(centroids = cents), 
                          type = c("centroids"),
                          molecules = data[[mol.type]], 
                          assay = assay)
      coords <- subset(x = coords, 
                       cells = intersect(x = Cells(x = coords[["centroids"]]),
                                         y = Cells(x = obj))) 
    }
  } else if ("segmentations" %in% names(data)) {
    segs <- CreateSegmentation(data[["segmentations"]])
    cents <- CreateCentroids(data[["centroids"]])
    segmentations.data <- list(centroids = cents, segmentation = segs)
    if (verbose) { 
      message("Creating FOVs..", "\n", 
              ">>> using segmentations") 
    }
    coords <- CreateFOV(coords = segmentations.data, 
                        type = c("segmentation", "centroids"), 
                        molecules = data[[mol.type]], 
                        assay = assay)
    # only consider the cells we have counts and a segmentation.
    # Cells which don't have a segmentation are probably found in other z slices.
    coords <- subset(x = coords,
                     cells = intersect(x = Cells(x = coords[["segmentation"]]),
                                       y = Cells(x = obj)))
  }
  
  # add z-stack index for cells
  if (add.zIndex) { obj$z <- data$zIndex %>% pull(z) }
  
  # add metadata vars
  if (verbose) { message(">>> adding metadata infos") }
  if (c("metadata" %in% names(data))) {
    metadata <- match.arg(arg = "metadata", choices = names(data), several.ok = TRUE)
    meta.vars <- names(data[[metadata]])
    for (i in meta.vars %>% seq) {
      obj %<>% AddMetaData(metadata = data[[metadata]][[meta.vars[i]]], 
                           col.name = meta.vars[i])
    }
  }
  
  # sanity on fov name
  fov %<>% gsub("_|-", ".", .)
  
  if (verbose) { message(">>> adding FOV") }
  obj[[fov]] <- coords
  
  ## filter - keep cells with counts > 0
  # helper function to return metadata
  callmeta <- function (object = NULL) { return(object@meta.data) }
  nCount <- grep("nCount", callmeta(obj) %>% names, value = TRUE)
  if (any(obj[[nCount]] == 0)) {
    if (verbose) { message(">>> filtering object - keeping cells with counts > 0") }
    obj %<>% subset(subset = !!base::as.symbol(nCount) > 0)
  } else { if (verbose) { message(">>> all counts are > 0") } }
  
  if (update.object) { 
    if (verbose) { message("Updating object:") 
      obj %<>% UpdateSeuratObject()
    } else { 
      obj %<>% 
        UpdateSeuratObject() %>% 
        suppressMessages() } }
  
  if (verbose) { message("Object is ready!") } 
  return(obj)
  
}