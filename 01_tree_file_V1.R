#generate the tree file
library("lidR")
library("raster")
library ("sf")
library("fasterize")
library("stringr")
library("dplyr")


wdir = "D:/scripts/data_prep_fgr/20240403_github/tree_location_v1/"
setwd(wdir)

# input tree species distribution layer - polygon shapefile
forest <- read_sf(dsn = wdir, layer = "forest")

# input layer for damaged areas for which height threshold is different
wind_areas <- read_sf(dsn = wdir, layer = "windthrown_areas")
#define the thresholds for tree height within and outside windthrown areas
thrs_wind <- 10; thrs_forest <- 5 

# input Canopy Height Model
CHM_f <- raster("CHM_filtered_input.tif"); plot(CHM_f)

# extent of the moving window for smoothing (median function) the chm
mw <- 3

# Set project CRS, check it (it has to be the same of dtm and treefile)
projectCRS = 32632L

#set allometric function to derive DBH from height values
DBH_fun <- function(height){ exp((height+22)/14) }   # Generic function


##### START OF THE ALGORITHM #####

#filtering extreme values >50 and <0 --> set them to 0
minValue(CHM_f);maxValue(CHM_f)
CHM_f <- reclassify(CHM_f, c(-Inf, 0, 0, 50, Inf, 50)); plot(CHM_f)
#writeRaster(CHM_f, "CHM_filtered.tif")


# create tile scheme
tilegrid = st_as_sf(terra::as.polygons(terra::ext(CHM_f))); plot(tilegrid, add=T)
tilegrid = st_as_sf(st_make_grid(tilegrid, cellsize = 1000)); plot(tilegrid, add=T)
raster.list = list()

i = 1 #for testing nrow(tilegrid)

y <- 1
for(i in 1:nrow(tilegrid) ){
  
  ss_orig = tilegrid[i,]
  ss = tilegrid[i,] %>%        
    st_buffer( 20, joinStyle="ROUND", endCapStyle = "FLAT")
  #plot(ss, add=TRUE)
  
  r.crop = crop(CHM_f, ss)[[1]] # crop raster
  #plot(r.crop)
  #r.crop <- projection(projectCRS)
  
  if(i == 1){ # add graphical progress
    plot(CHM_f)
    plot(tilegrid, add=TRUE) }
  
  plot(ss_orig, add=TRUE, col="light grey")
  
  
  # TEST TO NOT PROCESS TILES WITH ONLY NA VALUES
  val <- getValues(r.crop)
  cell_noNA <- sum(!is.na(val))
  
  if (cell_noNA>1000) {
    
    #START TREES IDENTIFICATION
    
    ### Smooth the CHM to better extract the tree tops and the crown radius
    CHM_f_smooth <- focal(as(r.crop, "SpatRaster"), w=mw, fun=median)
    CHM_f_smooth <- as(CHM_f_smooth, "Raster")
    
    ### START treetops and crown extraction from CHM
    
    write.table(locate_trees(CHM_f_smooth, lmf(ws=5, hmin = 5, shape = c("circular")))$geometry,row.names = F, col.names = F, "tops.txt")
    tree_tops <- read.delim("tops.txt", header = F)
    unlink("tops.txt")
    tree_tops <- as.data.frame(str_split_fixed(tree_tops$V1, "\\(|\\)|,", n = 5)); tree_tops$V1 <- NULL; tree_tops$V5 <- NULL
    colnames(tree_tops) <- c("x","y","z")
    options(digits=9); tree_tops$x <- as.numeric(tree_tops$x); tree_tops$y <- as.numeric(tree_tops$y); tree_tops$z <- as.numeric(tree_tops$z)
    ttops_points <- SpatialPointsDataFrame(coords = tree_tops, data = tree_tops)
    ttops_points$treeID = as.numeric(seq(1, nrow(ttops_points)))
    proj4string(ttops_points) <- crs(r.crop)
    
    crowns = silva2016(chm = CHM_f_smooth, treetops = ttops_points, exclusion = 0.3, max_cr_factor = 0.5)()
    contour = sf::as_Spatial(sf::st_as_sf(stars::st_as_stars(crowns), as_points = FALSE, merge = TRUE))
    contour$area <- area(contour)
    colnames(contour@data) <- c("treeID","area")
    cont = contour@data %>%  group_by(treeID) %>%  summarise(area = max(area))
    #write.table(cont, "cont.txt")
    #writeOGR(contour, dsn=wdir, layer="crowns_contour", driver = "ESRI Shapefile", overwrite_layer = T)
    trees_cr = merge(ttops_points, cont, by="treeID")
    trees_cr <- st_as_sf(trees_cr)
    
    ### END treetops and crown extraction from CHM
    
    
    # clipping tree tops in accordance with the original tile (not buffered)
    trees_cr <- st_transform(trees_cr,projectCRS)
    ss_orig <- st_set_crs(ss_orig,projectCRS)
    trees_cr <- st_intersection(trees_cr,ss_orig)
    trees_mat <- as.data.frame(trees_cr)
    
    
    ### START assigning tree species
    #assign to the tree points the species - species name is the same that will be passed to forestGALES
    forest_rast <- fasterize(forest, CHM_f_smooth, field="cat_num")
    spec <- extract(forest_rast,trees_cr)
    #plot(forest_rast); plot(trees_cr, add=T)
    trees_mat$spec <- spec
    
    trees_mat$spec_id <- "NS" #NA data -> spruce
    trees_mat$spec_id[trees_mat$spec==1] <- "NS" # silver fir -> 
    trees_mat$spec_id[trees_mat$spec==2] <- "BE" # Aceri-frassineto e aceri tiglieti
    trees_mat$spec_id[trees_mat$spec==3] <- "OK" # Alnete
    trees_mat$spec_id[trees_mat$spec==4] <- "NS" # Shrubs
    trees_mat$spec_id[trees_mat$spec==5] <- "BI" # Birch
    trees_mat$spec_id[trees_mat$spec==6] <- "BE" # Chestnut
    trees_mat$spec_id[trees_mat$spec==7] <- "BE" # Beech
    trees_mat$spec_id[trees_mat$spec==8] <- "NS" # Antropogenic forests
    trees_mat$spec_id[trees_mat$spec==9] <- "EL" # Larch and Larch with pinus cembrae
    trees_mat$spec_id[trees_mat$spec==10] <- "NS" # pinus mugo
    trees_mat$spec_id[trees_mat$spec==11] <- "OK" # Orno-ostrieti e ostrio-querceti
    trees_mat$spec_id[trees_mat$spec==12] <- "NS" # Spruce
    trees_mat$spec_id[trees_mat$spec==13] <- "NS" # Spruce + beech
    trees_mat$spec_id[trees_mat$spec==14] <- "SP" # Scots pine
    trees_mat$spec_id[trees_mat$spec==15] <- "OK" # Querco-carpineti e carpineti
    trees_mat$spec_id[trees_mat$spec==16] <- "NS" # Saliceti e altre formazioni riparie
    
    trees_mat$spec <- NULL
    
    treefile = trees_mat
    treefile$treeID <- NULL; treefile$x.1 <- NULL; treefile$y.1 <- NULL; treefile$z.1 <- NULL; treefile$geometry <- NULL
    colnames(treefile) = c("X","Y","H","area","species")
    ### END assigning tree species
    
    # 
    # #START different threshold height for defining a tree within and outside windthrown areas
    # trees_clean_int <- lengths(st_intersects(trees_cr,wind_areas))>0
    # trees_mat$wind <- trees_clean_int
    # trees_mat_thr <- subset(trees_mat, (wind==TRUE & z>=thrs_wind)|(wind==FALSE & z>=thrs_forest))
    # trees_mat_thr$z <- round(trees_mat_thr$z, digits=2)
    # trees_mat_thr$geometry <- NULL
    # #write.table(trees_mat_thr, "trees_ID_XYZ_cr_sp_wind.txt",col.names = T, row.names = F, sep="\t")
    # ttops_points <- as(st_as_sf(x = trees_mat_thr, coords = c("x", "y")), "Spatial") 
    # crs <- crs(CHM_f_smooth); crs(ttops_points)<- crs
    # #writeOGR(ttops_points, dsn=wdir, layer="trees_ID_XYZ_cr_sp_wind", driver = "ESRI Shapefile", overwrite_layer = T)
    # 
    # treefile = trees_mat_thr
    # treefile$treeID <- NULL; treefile$x.1 <- NULL; treefile$y.1 <- NULL; treefile$z.1 <- NULL;
    # colnames(treefile) = c("X","Y","H","area","species","wind")
    # #END different threshold height for defining a tree within and outside windthrown areas
    
    #START computing DBH from height
    # add DBH (calculated from the H/D curve)
    treefile$DBH = round(DBH_fun(treefile$H), digits = 2)
    
    #randomize calculation of DBH
    rand_DBH <- 0.10  # 10% of randomization
    num_rand <- runif(nrow(treefile), -rand_DBH, rand_DBH)
    treefile$DBH <- treefile$DBH + (treefile$DBH * num_rand)
    treefile$DBH <- round(treefile$DBH, digits = 2)
    #END computing DBH from height
    
    treefile$H = round(treefile$H, digits = 2)
    treefile$treeID = seq(1, nrow(treefile))
    treefile$date = format(Sys.Date(), "%d/%m/%Y")
    
    
    
    #Aggregating the treefiles of the different tiles
    if (y==1) {
      treefiletot <- treefile
    } else {
      treefiletot <- rbind(treefiletot,treefile)
    }
    
    y <- y+1
    
  }
  
  plot(ss_orig, add=TRUE, col="green") # positive progress
  
}

treefiletot$treeID = seq(1, nrow(treefiletot))

write.table(treefiletot, "trees_ID_XYZ_cr_sp_DBH.txt",col.names = T, row.names = F, sep="\t")

treefile_points <-  treefiletot
coordinates(treefile_points)=~X+Y
proj4string(treefile_points)<- crs(r.crop)
treefile_points <- st_as_sf(treefile_points)

st_write(treefile_points, dsn=wdir, layer="treefile_points", driver = "ESRI Shapefile", overwrite_layer = T)




