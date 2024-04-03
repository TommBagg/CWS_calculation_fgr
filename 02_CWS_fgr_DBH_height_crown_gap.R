### WIND RISK MAPPING - MUNICIPALITY LEVEL
### Authors: Tommaso Baggio, Maximiliano Costa, Niccolò Marchi

### Study area: Cordevole
### Scenario: Test pre VAIA

library("terra")
library("raster")
library("dplyr")
library("sf")
library("fgr")
library("fasterize")
library("OasisR")


wdir = "D:/scripts/data_prep_fgr/20240403_github/CWS_sim_v1/"
setwd(wdir)
out.folder = wdir

#Output raster name
out_name = "test_CWS_3X3.tif"

# Set the resolution of rasters, our idea is to set it at least to 20m
rast.resol = 20L

# Parameters for forest gap identification; minimum area for gap identification / minimum mean width for gap identification
min.gap.area = 1200; min.gap.width = 40

# Set the DLF threshold - if set to 1 the DLF is not considered - enable it to simulate the presence of snow
DLF.thr = 1

# snow height in m on the crown 
snow.height = 0.0

# choose the moving window size expressed in pixel between "3x3" and "5x5" for calculation of tree spacing and dominant height
MW.Spacing.Hmax = "3x3"

# Set project CRS, check it (it has to be the same of dtm and treefile)
projectCRS = 32632L # for Rocca Pietore we are using WGS84/UTM zone 32N, EPSG:32632

# import species parameters (customized)
load(file = "mydata1.sp.rda")

# import fgr constant (customized, DLF = 1)
load(file = "mydata.fc.rda")

# import raster for area of interest and to derive the slope (classification of NS parameters)
r = rast(file.path("D:/progetti/resilience/Gis/dtm_5m/DTM_5m_cl.tif")); plot(r)
crs(r) <- paste0("epsg:",projectCRS)

# Import file of tree locations and filter out based on tree characteristics
treefile <- read.delim("D:/progetti/resilience/fgr_cordevole/00_tree_tops_alg_v1/trees_ID_XYZ_cr_sp_DBH.txt")
## Max tree heigth supposed equal to 40m ######################
treefile <- subset(treefile, H<40)
## Minimum value for the crown area set to > 3m^2 (equivalent to a diameter of 1.95) ######################
treefile <- subset(treefile, area>4)
## Use the improved species specific parameters for european larch
treefile$species[treefile$species=="EL"] <- "U_EL_t3"

# convert to SPDF keeping X & Y 
treefile = st_as_sf(treefile, coords = c("X","Y"))
treefile = subset(treefile, !is.na(DBH)) #  remove rows for which DBH==NA
#st_write(treefile, "FGsim_tree_tops.shp", driver="ESRI Shapefile")



#++++++++++++++++++++++++++++++++++++++++++++++++++++

#### Rasterise function ####

## for testing input data
#treefile.temp = trees.crop; r.temp = r.crop; method = "tmc"; resol=20L; standID.field="standID"; H.field="H"; DBH.field="DBH"; area.field="area"; MW_Spacing_Hmax="3x3"; min_gap_area=1200; min_gap_width=40; DLF_thr=1; snow_height=0
##

rasterise_fgr = function(treefile.temp, r.temp, forest, standID.field="standID", H.field="H", DBH.field="DBH", area.field="area",
                         method=c("tmc","roughness"), resol=20, snow_height=0, DLF_thr=1, min_gap_area=1200, min_gap_width=40, MW_Spacing_Hmax="3x3"){
  
  # safety checks
  if( !is(treefile.temp,"SpatialPointsDataFrame") & !is(treefile.temp,"sf") ){
    stop("Treefile is not a spatial object")
  } else { if( is(treefile.temp,"SpatialPointsDataFrame") ){ treefile.temp = st_as_sf(treefile.temp) } }
  if( is.na(match(H.field, names(treefile.temp))) | is.na(match(DBH.field, names(treefile.temp))) ){ stop("One or more fields were not found: check spelling")}
  if( is.na(crs(r.temp)) ){ stop("raster's CRS is missing")}
  if( resol < 20 ){ warning("Resolution is too small: results might loose meaning")}
  stopifnot(is.character(method))
  # add a check for CRS
  
  r.temp = terra::terrain(r.temp, v="slope", unit="degrees") # calculate slope
  treefile.temp$slope = round(terra::extract(r.temp, vect(treefile.temp))[,2]) # assign slope to trees
  r.temp2 <- r.temp
  
  # remove points outside study area (-> with slope == NA)
  treefile.temp = treefile.temp[!is.na(treefile.temp$slope),]
  #warning("Trees outside the raster extent have been removed")
  
  if( nrow(trees.crop)==0 ){
    
    out = r.temp
    out[out] = NA
    
    out = c(out,out,out) # stack layers
    return(out)
    
  } else {
    
    ### Select species parameters according to dbh or slope, based on fgr database and italian pulling tests
    
    treefile.temp$species[treefile.temp$species=="NS" & treefile.temp$DBH > 40 & treefile.temp$slope > 20] <- "U_NS_s"
    treefile.temp$species[treefile.temp$species=="NS" & treefile.temp$DBH > 40 & treefile.temp$slope <= 20] <- "U_NS_f"
    treefile.temp$species[treefile.temp$species=="NS" & treefile.temp$DBH <= 40 & treefile.temp$slope > 20] <- "U_NS_fgr_s"
    treefile.temp$species[treefile.temp$species=="NS" & treefile.temp$DBH <= 40 & treefile.temp$slope <= 20] <- "NS"
    
    if( floor(resol/res(r.temp)[1]) != 1 ){ # aggregate only if resolutions are different
      r.temp = terra::aggregate(r.temp, fact=floor(resol/res(r.temp)[1]), fun=mean, na.rm=TRUE)# resample to coarser resolution
    }
    
    # create a pixel ID from coords
    #r.temp = setValues(r.temp, factor(paste(crds(r.temp)[,1]),round(crds(r.temp)[,2]),sep="_"))) # there are issues with factors in raster
    #r.temp1 = set.values(r.temp, cells=seq(1,ncell(r.temp),1), values=seq(1,ncell(r.temp),1) ) # simply use a progressive number
    r.temp = setValues(r.temp, seq(1,ncell(r.temp),1))
    
    treefile.temp$treeID = terra::extract(r.temp, vect(treefile.temp))[,2] # assign pixel ID to trees
    
    r.df = as.data.frame(r.temp, xy=TRUE)[,1:3]
    colnames(r.df) = c("X","Y","treeID")
    
    #create a copy for Hmax calculation
    r.temp1 <- r.temp
    treefile.temp1 <- treefile.temp
    
    values(r.temp) = NA # make sure the raster template is empty
    
    # summarise and assign stand-level values: Hmean, Hmax
    # https://stackoverflow.com/questions/29678435/how-to-pass-dynamic-column-names-in-dplyr-into-custom-function
    
    tmp = treefile.temp %>% # species prevalence
      st_drop_geometry() %>%
      group_by(treeID) %>%
      count(species) %>%
      slice_max(n,n=1) %>%
      dplyr::select( -n)
    
    treefile.temp = treefile.temp %>%
      group_by(treeID) %>%
      summarise(px_Hmean = mean(!!as.name(H.field)),
                px_DBHmean = mean(!!as.name(DBH.field)),
                px_CRWIDTHmean = mean(!!as.name(area.field)))

    treefile.temp$px_CRWIDTHmean <- ((treefile.temp$px_CRWIDTHmean/pi)^0.5)*2
    treefile.temp = merge(treefile.temp, tmp, by="treeID")

    
    
    ### MOVING WINDOW 3x3
    if (MW_Spacing_Hmax=="3x3") {

      ### START compute pxHmax with a moving window of 3x3, 90th percentile ###
      treefile.temp$px_Hmax <- 0.0
      treefile.temp <- treefile.temp[order(treefile.temp$treeID),]
      treefile.temp1 <- treefile.temp1[order(treefile.temp1$treeID),]

      trees_h <- data.frame(treefile.temp1$H,treefile.temp1$treeID)
      colnames(trees_h) <- c("H","id")
      trees_h <- trees_h[order(trees_h$id),]
      ncol_r.temp1 <- ncol(r.temp1)

      j <-1
      for (i in c(treefile.temp$treeID)) {
        N <- subset(trees_h, trees_h$id==i-ncol_r.temp1+1 | trees_h$id==i-ncol_r.temp1 | trees_h$id==i-ncol_r.temp1-1)[,1]
        C <- subset(trees_h, trees_h$id==i-1 | trees_h$id==i | trees_h$id==i+1)[,1]
        S <- subset(trees_h, trees_h$id==i+ncol_r.temp1+1 | trees_h$id==i+ncol_r.temp1 | trees_h$id==i+ncol_r.temp1-1)[,1]
        H3x3 <- c(N,C,S)
        treefile.temp$px_Hmax[j] <- quantile(H3x3, 0.90)
        j <- j+1
      }

      ### END compute pxHmax


      ### START compute px_spacing with a moving window of 3x3
      treefile.temp$px_spacing <- 0.0
      trees_id <- trees_h$id

      #i<-35
      j <-1
      for (i in c(treefile.temp$treeID)) {
        N <- NA; NE <- NA; E <- NA; SE <- NA; S <- NA; SW <- NA; W <- NA; NW <- NA; C <- NA

        N <- length(trees_id[trees_id==i-ncol_r.temp1])
        NE <- length(trees_id[trees_id==i-ncol_r.temp1+1])
        NW <- length(trees_id[trees_id==i-ncol_r.temp1-1])
        S <- length(trees_id[trees_id==i+ncol_r.temp1])
        SE <- length(trees_id[trees_id==i+ncol_r.temp1+1])
        SW <- length(trees_id[trees_id==i+ncol_r.temp1-1])
        C <- length(trees_id[trees_id==i])
        E <- length(trees_id[trees_id==i+1])
        W <- length(trees_id[trees_id==i-1])

        num <- N+NE+NW+S+SE+SW+C+W+E

        if (N>0) {N<-1}
        if (NE>0) {NE<-1}
        if (NW>0) {NW<-1}
        if (S>0) {S<-1}
        if (SE>0) {SE<-1}
        if (SW>0) {SW<-1}
        if (C>0) {C<-1}
        if (E>0) {E<-1}
        if (W>0) {W<-1}
        num_cell <- N+NE+NW+S+SE+SW+C+W+E
        treefile.temp$px_spacing[j] <- sqrt((res(r.temp)[1]^2)*num_cell/num)
        j <- j+1

      }

      ### END compute px_spacing

    }

    ### MOVING WINDOW 5x5
    if (MW_Spacing_Hmax=="5x5") {

      ### compute pxHmax with a moving window of 5x5, 90th percentile
      treefile.temp$px_Hmax <- 0.0
      treefile.temp <- treefile.temp[order(treefile.temp$treeID),]
      treefile.temp1 <- treefile.temp1[order(treefile.temp1$treeID),]

      trees_h <- data.frame(treefile.temp1$H,treefile.temp1$treeID)
      colnames(trees_h) <- c("H","id")
      trees_h <- trees_h[order(trees_h$id),]
      ncol_r.temp1 <- ncol(r.temp1)

      j <-1
      for (i in c(treefile.temp$treeID)) {
        NN <- subset(trees_h, trees_h$id==i-2*ncol_r.temp1+2 | trees_h$id==i-2*ncol_r.temp1+1 | trees_h$id==i-ncol_r.temp1*2 | trees_h$id==i-ncol_r.temp1*2-1 | trees_h$id==i-ncol_r.temp1*2-2)[,1]
        N <- subset(trees_h, trees_h$id==i-ncol_r.temp1+2 | trees_h$id==i-ncol_r.temp1+1 | trees_h$id==i-ncol_r.temp1 | trees_h$id==i-ncol_r.temp1-1 | trees_h$id==i-ncol_r.temp1-2)[,1]
        C <- subset(trees_h, trees_h$id==i-2 | trees_h$id==i-1 | trees_h$id==i | trees_h$id==i+1 | trees_h$id==i+2)[,1]
        S <- subset(trees_h, trees_h$id==i+ncol_r.temp1+2 | trees_h$id==i+ncol_r.temp1+1 | trees_h$id==i+ncol_r.temp1 | trees_h$id==i+ncol_r.temp1-1 | trees_h$id==i+ncol_r.temp1-2)[,1]
        SS <- subset(trees_h, trees_h$id==i+2*ncol_r.temp1+2 | trees_h$id==i+2*ncol_r.temp1+1 | trees_h$id==i+ncol_r.temp1*2 | trees_h$id==i+ncol_r.temp1*2-1 | trees_h$id==i+ncol_r.temp1*2-2)[,1]
        H5x5 <- c(NN,N,C,S,SS)
        treefile.temp$px_Hmax[j] <- quantile(H5x5, 0.90)
        j <- j+1
      }
      ### end compute pxHmax


      ### compute px_spacing with a moving window of 5x5
      treefile.temp$px_spacing <- 0.0
      trees_id <- trees_h$id

      # A B   C    D E    #MW 5X5
      # F G   H    I L
      # M N CENTER O P
      # Q R   S    T U
      # V W   X    Y Z

      #i<-35
      j <-1
      for (i in c(treefile.temp$treeID)) {

        A <- NA; B <- NA; C <- NA; D <- NA; E <- NA
        FF <- NA; G <- NA; H <- NA; I <- NA; L <- NA
        M <- NA; N <- NA; CENTTER <- NA; O <- NA; P <- NA
        Q <- NA; R <- NA; S <- NA; TT <- NA; U <- NA
        V <- NA; W <- NA; X <- NA; Y <- NA; Z <- NA

        A <- length(trees_id[trees_id==i-ncol_r.temp1*2-2])
        B <- length(trees_id[trees_id==i-ncol_r.temp1*2-1])
        C <- length(trees_id[trees_id==i-ncol_r.temp1*2])
        D <- length(trees_id[trees_id==i-ncol_r.temp1*2+1])
        E <- length(trees_id[trees_id==i-ncol_r.temp1*2+2])

        FF <- length(trees_id[trees_id==i-ncol_r.temp1-2])
        G <- length(trees_id[trees_id==i-ncol_r.temp1-1])
        H <- length(trees_id[trees_id==i-ncol_r.temp1])
        I <- length(trees_id[trees_id==i-ncol_r.temp1+1])
        L <- length(trees_id[trees_id==i-ncol_r.temp1+2])

        M <- length(trees_id[trees_id==i-2])
        N <- length(trees_id[trees_id==i-1])
        CENTER <- length(trees_id[trees_id==i])
        O <- length(trees_id[trees_id==i+1])
        P <- length(trees_id[trees_id==i+2])

        Q <- length(trees_id[trees_id==i+ncol_r.temp1-2])
        R <- length(trees_id[trees_id==i+ncol_r.temp1-1])
        S <- length(trees_id[trees_id==i+ncol_r.temp1])
        TT <- length(trees_id[trees_id==i+ncol_r.temp1+1])
        U <- length(trees_id[trees_id==i+ncol_r.temp1+2])

        V <- length(trees_id[trees_id==i+ncol_r.temp1*2-2])
        W <- length(trees_id[trees_id==i+ncol_r.temp1*2-1])
        X <- length(trees_id[trees_id==i+ncol_r.temp1*2])
        Y <- length(trees_id[trees_id==i+ncol_r.temp1*2+1])
        Z <- length(trees_id[trees_id==i+ncol_r.temp1*2+2])

        num <- A+B+C+D+E+FF+G+H+I+L+M+N+CENTER+O+P+Q+R+S+TT+U+V+W+X+Y+Z

        if (A>0) {A<-1}; if (FF>0) {FF<-1}; if (M>0) {M<-1}; if (Q>0) {Q<-1}; if (V>0) {V<-1}
        if (B>0) {B<-1}; if (G>0) {G<-1}; if (N>0) {N<-1}; if (R>0) {R<-1}; if (W>0) {W<-1}
        if (C>0) {C<-1}; if (H>0) {H<-1}; if (CENTER>0) {CENTER<-1}; if (S>0) {S<-1}; if (X>0) {X<-1}
        if (D>0) {D<-1}; if (I>0) {I<-1}; if (O>0) {O<-1}; if (TT>0) {TT<-1}; if (Y>0) {Y<-1}
        if (E>0) {E<-1}; if (L>0) {L<-1}; if (P>0) {P<-1}; if (U>0) {U<-1}; if (Z>0) {Z<-1}

        num_cell <- A+B+C+D+E+FF+G+H+I+L+M+N+CENTER+O+P+Q+R+S+TT+U+V+W+X+Y+Z
        treefile.temp$px_spacing[j] <- sqrt((res(r.temp)[1]^2)*num_cell/num)
        j <- j+1

      }

      ### END compute px_spacing

    }
    
    
    # ### START distance to the edge and gap size
    # p <- as_Spatial(sf::st_as_sf(as.points(r.temp1))); colnames(p@data) <- "treeID"
    # treefile.temp <- merge(treefile.temp,p , by="treeID")
    # treefile.temp$geometry <- NULL; 
    # treefile.temp = treefile.temp %>% rename(x = coords.x1)
    # treefile.temp = treefile.temp %>% rename(y = coords.x2)
    # coordinates(treefile.temp) <- ~x + y
    # treefil <- as_Spatial(sf::st_as_sf(treefile.temp))
    # r.trees = terra::rasterize(vect(treefil), r.temp,background=0)#; plot(r.trees)
    # r.trees[r.trees==1] <- NA
    # r.trees_buf_v <- as.polygons(r.trees, dissolve = T)#; plot(r.trees_buf_v, add=T)
    # #polya <- as(sf::st_as_sf(r.trees_buf_v), "Spatial")
    # polybd <- terra::disagg((r.trees_buf_v))
    # 
    # # calculation of the width and length of the gap features
    # 
    # polybd$area <- area(polybd)#; plot(polybd, col="red")
    # polybd$perimeter <- perimeter(polybd)
    # polybd$width1 <- (polybd$perimeter - (polybd$perimeter^2 - 16 * polybd$area)^0.5)/4
    # polybd$width2 <- (polybd$perimeter + (polybd$perimeter^2 - 16 * polybd$area)^0.5)/4
    # polybd$size_fin <- (polybd$width1 + polybd$width2)/2
    # # here improve the polygon representing the forest - try to exclude feature within the forest
    # 
    # polybd <- as(sf::st_as_sf(polybd), "Spatial")
    # polybd.sub <- subset(polybd, area > min_gap_area) # removed features with area lower than threshold
    # polybd.sub <- subset(polybd.sub, width1 > min_gap_width)#; plot(polybd.sub, col="red") # removed features with width lower than 40 m (around 2*mean tree height)
    # 
    # treefil <- vect(treefil); crs(treefil) <- paste0("epsg:",projectCRS)
    # polybd.sub <- vect(polybd.sub); crs(polybd.sub) <- paste0("epsg:",projectCRS)
    # 
    # #extract forest points and derive point polygon (representing non-forest areas) distance
    # mat_dist <- terra::distance(treefil, polybd.sub, unit="m", pairwise=T)
    treefile.temp$px_dist_edge <- NA
    treefile.temp$gap_size <- NA
    # plot(polybd.sub); plot(treefil, add=T)
    # # i <- 1
    # for (i in 1:ncol(mat_dist)) {
    #   treefile.temp$px_dist_edge[i] <- min(mat_dist[,i])
    #   row_pol <- match(treefile.temp$px_dist_edge[i],mat_dist[,i])
    #   treefile.temp$gap_size[i] <- polybd.sub$size_fin[row_pol]
    # }
    #View(treefile.temp@data)
    #writeVector(vect(polybd.sub), "test2")
    #writeRaster(out, "out4.tif")
    
    ### END distance to the edge and gap size
    
    ### START assign different DLF thr in accordance to forest species
    
    treefile.temp$DLF_thr <- DLF_thr
    #for norway spruce assign to 1 because calibrated by tree pulling experiments performed at DBH
    #treefile.temp$DLF_thr[treefile.temp$species=="U_NS_s"] <- 1
    #treefile.temp$DLF_thr[treefile.temp$species=="U_NS_f"] <- 1
    #treefile.temp$DLF_thr[treefile.temp$species=="U_NS_fgr_s"] <- 1
    #treefile.temp$DLF_thr[treefile.temp$species=="NS"] <- 2.5
    
    ### END assign DLF thr
    
    # assign SNOW HEIGHT
    treefile.temp$snow_height = snow_height
    # END assign snow height
    
    treefile.temp$standID = "1" # add fake standID
    
    if( "tmc" %in% method ){
      
      # initialise the model
      args.to.vectorize <- formalArgs(fg_tmc)
      args.to.vectorize <- args.to.vectorize[!args.to.vectorize %in% c("species_parameters", "fgr_constants")]
      
      # calculate using the TMC Method
      fg_tmc_v <- Vectorize(fg_tmc, vectorize.args = args.to.vectorize)
      
      AA.stand_id <- as.character(treefile.temp$standID)
      AA.tree_id <- as.character(treefile.temp$treeID)
      AA.date <- 'X'
      AA.species <- as.character(treefile.temp$species)
      AA.tree_ht <- treefile.temp$px_Hmean
      AA.dbh <- treefile.temp$px_DBHmean
      AA.spacing_current <- treefile.temp$px_spacing
      AA.stand_mean_ht <- treefile.temp$px_Hmean
      AA.stand_top_ht <- treefile.temp$px_Hmax
      AA.stand_mean_dbh <- treefile.temp$px_DBHmean
      AA.cr_width <- treefile.temp$px_CRWIDTHmean
      AA.stand_cr_width <- treefile.temp$px_CRWIDTHmean
      AA.dist_edge <- treefile.temp$px_dist_edge
      AA.predominant_species <- as.character(treefile.temp$species)
      AA.gap_size <- treefile.temp$gap_size
      AA.dlf_thr <- treefile.temp$DLF_thr
      AA.snow_height =treefile.temp$snow_height
      AA.full_output <- 0
      
      
      out = as.data.frame(t(fg_tmc_v(stand_id = AA.stand_id, tree_id = AA.tree_id, date = AA.date,
                                     species = AA.species, tree_ht = AA.tree_ht, dbh = AA.dbh,
                                     spacing_current = AA.spacing_current, predominant_species = AA.predominant_species,
                                     stand_mean_dbh = AA.stand_mean_dbh, stand_top_ht = AA.stand_top_ht,
                                     cr_width = AA.cr_width, stand_cr_width = AA.stand_cr_width, dist_edge = AA.dist_edge,
                                     full_output = AA.full_output, species_parameters = mydata.sp, gap_size = AA.gap_size,
                                     dlf_threshold = AA.dlf_thr, snow_depth = AA.snow_height, fgr_constants = mydata.fc)))

      
      
      
      out = tidyr::unnest(out, cols =  c(stand_id, tree_id, date, species, u_elev_damage, mode_of_damage,
                                         u_elev_b, u_elev_o, prob_damage, prob_b, prob_o, Warnings))
      
      names(out)[grep("tree", names(out))] = 'treeID'
      
    } else if( "roughness" %in% method ){
      
      # initialize the model
      args.to.vectorize <- formalArgs(fg_rou)
      args.to.vectorize <- args.to.vectorize[!args.to.vectorize %in% c("species_parameters", "fgr_constants")]
      
      # calculate using the Roughness Method
      args.to.vectorize <- formalArgs(fg_rou)
      args.to.vectorize <- args.to.vectorize[!args.to.vectorize %in% c("species_parameters", "fgr_constants")]
      fg_rou_v <- Vectorize(fg_rou, vectorize.args = args.to.vectorize)
      
      AA.stand_id <- as.character(treefile.temp$treeID)
      AA.date <- 'X'
      AA.species <- as.character(treefile.temp$species)
      AA.mean_ht <- treefile.temp$px_Hmean
      AA.mean_dbh <- treefile.temp$px_DBHmean
      AA.spacing <- treefile.temp$px_spacing
      AA.full_output <- 0
      
      out = as.data.frame(t(fg_rou_v(AA.stand_id, AA.date, AA.species, AA.mean_ht, AA.mean_dbh,
                                     AA.spacing, AA.full_output, species_parameters = mydata.sp,
                                     fgr_constants = mydata.fc)))
      
      out <- tidyr::unnest(out, cols = c(stand_id, date, species, u_elev_damage, mode_of_damage, u_elev_b,
                                         u_elev_o, prob_damage, prob_b, prob_o, Warning_Variables,
                                         Max_Stem_Weight_Warning))
      
      names(out)[grep("stand", names(out))] = 'treeID'
      
    } else { stop("Wrong method passed")}
    
    out = merge(r.df, out, by="treeID")
    
    out = st_as_sf(out, coords = c("X","Y")#,
                   #proj4string = crs(r.temp)
    )
    
    # rasterise variables
    r1 = rasterize(vect(out), r.temp, "u_elev_damage") #
    
    r2 = rasterize(vect(out), r.temp, "u_elev_b") #
    
    r3 = rasterize(vect(out), r.temp, "u_elev_o") #
    
    # output a stack containing the two varbs
    out = c(r1, r2, r3)
    names(out) = c("u_elev_damage","u_elev_b","u_elev_o")
    
    return( out )
    
  }
  
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#### ANALYSIS ####

# create tile scheme
tilegrid = st_as_sf(terra::as.polygons(terra::ext(r)))
tilegrid = st_as_sf(st_make_grid(tilegrid, cellsize = 1000)); plot(tilegrid, add=T)
raster.list = list()

tictoc::tic()

#### To solve "Error: [`[<-`] lengths of cells and values do not match" try changing the buffer value of the tile
i = 1
for(i in 1:nrow(tilegrid) ){
  
  ss_orig = tilegrid[i,]
  ss = tilegrid[i,] %>%       # buffer length increased to calculate the length of the upwind edge (suggested to 10*mean_height) 
    st_buffer( 50, joinStyle="ROUND", endCapStyle = "FLAT") # da correggere l'angolo e anche valore del buffer se da errore
  
  r.crop = crop(r, ss)[[1]] # crop raster
  trees.crop = st_intersection(treefile, ss) # crop treetops
  
  if(i == 1){ # add graphical progress
    plot(r)
    plot(tilegrid, add=TRUE) }
  
  plot(ss_orig, add=TRUE, col="light grey")
  
  if( nrow(trees.crop)==0 ){
    r.final = r.crop
    r.final = terra::aggregate(r.final, fact=floor(rast.resol/res(r.final)[1]), fun=mean, na.rm=TRUE)# resample to coarser resolution
    values(r.final) = NA
  } else {
    r.final = rasterise_fgr(treefile.temp = trees.crop , r.temp = r.crop, method = "tmc", resol = rast.resol, snow_height=snow.height, DLF_thr=DLF.thr, min_gap_area=min.gap.area, min_gap_width=min.gap.width, MW_Spacing_Hmax=MW.Spacing.Hmax) # run the Rasterise function
    #r.final = rasterise_fgr(trees.crop, r.temp = r.crop, method = "roughness") # run the Rasterise function
  }
  r.final = crop(r.final, ss_orig)
  raster.list[[paste0("tile_",i)]] = r.final # add to list
  plot(ss_orig, add=TRUE, col="green") # positive progress
}

tictoc::toc()

final.raster = terra::mosaic(terra::sprc(raster.list))

plot(final.raster[[1]])
plot(final.raster[[2]])
plot(final.raster[[3]])

#Check the name of the output file
writeRaster(final.raster, filename=file.path(out.folder, out_name), overwrite=TRUE, wopt=list(gdal=c("COMPRESS=LZW", "TFW=YES","of=GTiff")))





