Cluster.plot <- function(method, parameter, GIS.Age, country,GADM,GIS.level){
  switch(method,
         "moran" = {
           #get out island district
           moran_poly <- GADM[[3]][-c(125, 134, 145, 153, 158, 161, 162, 189, 190, 193, 195),]


           switch(GIS.Age,
                  "no"={
                    df_moran <- CDM.table[, c("crd_sir", "gadm_id", "outcome_count", "crd_expected")]
                    colnames(df_moran) <- c("SIR", "gadm_id", "outcome_count", "expected")

                  },
                  "yes"={
                    df_moran <- CDM.table[, c("std_sir", "gadm_id", "outcome_count", "std_expected")]
                    colnames(df_moran) <- c("SIR", "gadm_id", "outcome_count", "expected")

                  }

           )

           moran_poly@data <- moran_poly@data %>%
             left_join(dplyr::select(df_moran, SIR, gadm_id, outcome_count, expected), by=c("ID_2"="gadm_id"))

           spatmatrix <- poly2nb(moran_poly)

           # create a neighbours list with spatial weights
           listw <- nb2listw(spatmatrix)

           # calculate the local moran of the distribution of target population
           lmoran <- localmoran(moran_poly$SIR, listw)
           summary(lmoran)

           # padronize the variable and save it to a new column
           moran_poly$SIR_s <- scale(moran_poly$SIR)  %>% as.vector()

           # create a spatially lagged variable and save it to a new column
           moran_poly$lag_s_SIR <- lag.listw(listw, moran_poly$SIR_s)

           # summary of variables, to inform the analysis
           summary(moran_poly$SIR_s)
           summary(moran_poly$lag_s_SIR)

           # moran sccaterplot, in basic graphics (with identification of influential observations)
           x <- moran_poly$SIR_s
           y <- moran_poly$lag_s_SIR %>% as.vector()
           xx <- data_frame(x, y)

           moran_poly$quad_sig <- NA

           # high-high quadrant
           moran_poly[(moran_poly$SIR_s >= 0 &
                         moran_poly$lag_s_SIR >= 0) &
                        (lmoran[, 5] <= 0.05), "quad_sig"] <- "high-high"
           # low-low quadrant
           moran_poly[(moran_poly$SIR_s <= 0 &
                         moran_poly$lag_s_SIR <= 0) &
                        (lmoran[, 5] <= 0.05), "quad_sig"] <- "low-low"
           # high-low quadrant
           moran_poly[(moran_poly$SIR_s >= 0 &
                         moran_poly$lag_s_SIR <= 0) &
                        (lmoran[, 5] <= 0.05), "quad_sig"] <- "high-low"
           # low-high quadrant
           moran_poly@data[(moran_poly$SIR_s <= 0
                            & moran_poly$lag_s_SIR >= 0) &
                             (lmoran[, 5] <= 0.05), "quad_sig"] <- "low-high"
           # non-significant observations
           moran_poly@data[(lmoran[, 5] > 0.05), "quad_sig"] <- "not signif."

           moran_poly$quad_sig <- as.factor(moran_poly$quad_sig)
           moran_poly@data$id <- rownames(moran_poly@data)

           df <- fortify(moran_poly, region="id")
           df <- left_join(df, moran_poly@data)

           map <- GIS.background(GADM[[3]]@bbox, country)

           plot <- map +
             geom_polygon(data=df,aes(x=long,y=lat,group=group,fill=quad_sig),colour="black", size=.05)+
             coord_equal() +
             theme_void() + scale_fill_brewer( palette = "Set1") +
             coord_fixed(ratio=1.1)

         },
         "kulldorff" = {
           parameter <- as.numeric(parameter)

           gadm_id <- CDM.table[,"gadm_id"]

           ls <- list()
           for(i in 1:length(gadm_id)) {
             idx <- gadm_id[i]
             geo<-GADM[[3]]@polygons[[idx]]
             a<-idx
             y<-geo@labpt[1]
             x<-geo@labpt[2]
             theta <- cbind(a,x,y)
             theta <- data.frame(theta)
             ls[[i]] <- theta
           }


           df <- do.call(rbind, ls)

           geo <- SpatialEpi::latlong2grid(df[, c(2,3)])

           pop.upper.bound <- parameter
           n.simulations <- 999
           alpha.level <- 0.05
           n.strata <- 16
           plot <- TRUE

           df <- dplyr::left_join(CDM.table, df, by=c("gadm_id" = "a"))

           population <- tapply(df$target_count,df$gadm_id,sum)
           cases <- tapply(df$outcome_count,df$gadm_id,sum)

           switch(GIS.Age,
                  "no"={
                    df <- df[, c("gadm_id", "target_count", "outcome_count", "crd_expected", "x", "y")]
                    colnames(df) <- c("gadm_id", "target_count", "Observed", "Expected", "x", "y")
                    expected.cases <- df$Expected
                  },
                  "yes"={
                    df <- df[, c("gadm_id", "target_count", "outcome_count", "std_expected", "x", "y")]
                    colnames(df) <- c("gadm_id", "target_count", "Observed", "Expected", "x", "y")
                    expected.cases <- df$Expected
                  }
           )

           result <- SpatialEpi::kulldorff(geo, cases, population, expected.cases, pop.upper.bound,
                                           n.simulations, alpha.level, plot)
           n.cluster <- 1 + length(result$secondary.clusters)

           idxNum <- paste0("ID_", GIS.level)
           idxName <- paste0("NAME_", GIS.level)
           tempGADM <- GADM[[GIS.level+1]]
           tempGADM@data <- dplyr::left_join(GADM[[GIS.level+1]]@data, CDM.table, by = structure(names = idxNum ,"gadm_id"))

           cluster.indexs <- list(result$most.likely.cluster)
           if(!is.null(result$secondary.clusters)){
             for (i in 1:length(result$secondary.clusters)){
               cluster.indexs <- c(cluster.indexs,list(result$secondary.clusters[[i]]))
             }
           }


           tempGADM$k.cluster <- NA
           tempGADM$population <- NA
           tempGADM$number.of.cases <- NA
           tempGADM$expected.cases <- NA
           tempGADM$SMR <- NA
           tempGADM$log.likelihood.ratio <- NA
           tempGADM$monte.carlo.rank <- NA
           tempGADM$p.value <- NA
           for (i in 1:length(cluster.indexs)){
             temp <- cluster.indexs[[i]][[1]]
             tempGADM@data$k.cluster[temp] <- i
             tempGADM@data$population[temp] <- cluster.indexs[[i]]$population
             tempGADM@data$number.of.cases[temp] <- cluster.indexs[[i]]$number.of.cases
             tempGADM@data$expected.cases[temp] <- cluster.indexs[[i]]$expected.cases
             tempGADM@data$SMR[temp] <- cluster.indexs[[i]]$SMR
             tempGADM@data$log.likelihood.ratio[temp] <- cluster.indexs[[i]]$log.likelihood.ratio
             tempGADM@data$monte.carlo.rank[temp] <- cluster.indexs[[i]]$monte.carlo.rank
             tempGADM@data$p.value[temp] <- cluster.indexs[[i]]$p.value
           }


           polygon_popup <- paste0("<strong>Name: </strong>", tempGADM@data[, idxName], "<br>",
                                   "<strong>Target: </strong>", tempGADM@data$target_count, "<br>",
                                   "<strong>Outcome: </strong>", tempGADM@data$outcome_count, "<br>",
                                   "<strong>SIR: </strong>", round(tempGADM@data$std_sir, 2), " (", round(tempGADM@data$std_sirlower, 2), "-", round(tempGADM@data$std_sirupper, 2), ")", "<br>",
                                   "<strong>Proportion: </strong>", round(tempGADM@data$std_prop, 2), " (", round(tempGADM@data$std_proplower, 2), "-", round(tempGADM@data$std_propupper, 2), ")",
                                   "<br>------------------Clustering Result>------------------<br>",
                                   "<strong>population: </strong>", tempGADM@data$population,"<br>",
                                   "<strong>number.of.cases: </strong>", tempGADM@data$number.of.cases,"<br>",
                                   "<strong>expected.cases: </strong>", round(tempGADM@data$expected.cases,4),"<br>",
                                   "<strong>SMR: </strong>", round(tempGADM@data$SMR,4),"<br>",
                                   "<strong>log.likelihood.ratio: </strong>", round(tempGADM@data$log.likelihood.ratio,4),"<br>",
                                   "<strong>p.value: </strong>", round(tempGADM@data$p.value,3),"<br>"
           )

           max_value <- max(na.omit(unique(tempGADM@data$k.cluster)))
           pal <- colorNumeric(topo.colors(max_value),domain =  NULL, na.color = "#FFFFFF", alpha = FALSE,
                               reverse = FALSE)

           m <- leaflet(GADM[[GIS.level+1]]) %>%
             addTiles %>%
             fitBounds (
               lng1=GADM[[1]]@bbox[1,1], lng2=GADM[[1]]@bbox[1,2],
               lat1=GADM[[1]]@bbox[2,1], lat2=GADM[[1]]@bbox[2,2])

           #create leaflet map
           plot <- m %>% addPolygons(data = tempGADM,
                                     fillColor= pal(tempGADM@data$k.cluster),
                                     fillOpacity = 0.5,
                                     weight = 1,
                                     color = "black",
                                     dashArray = "3",
                                     popup = polygon_popup,
                                     highlight = highlightOptions(
                                       weight = 5,
                                       color = "#666",
                                       dashArray = "",
                                       fillOpacity = 0.7,
                                       bringToFront = TRUE)) %>%
             addLegend(pal = pal, values = ~tempGADM@data$k.cluster, opacity = 0.7, title = NULL,
                       position = "bottomright")

         }
  )
  return(plot)
}
