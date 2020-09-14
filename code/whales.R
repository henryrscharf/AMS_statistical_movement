## libraries ----
library(ggmap)
library(mgcv)
library(momentuHMM)
library(crawl)
library(ncdf4)
library(raster)
library(sf)
library(parallel)
## projection, colors, plot characteristics, maps, etc. ----
st_CRS <- CRS("+proj=utm +zone=11 +ellps=WGS84 +datum=WGS84 +units=km")
colors <- c(obs = "black", true = "black")
argos_colors <- c(rep("black", 4), rep("black", 2))
argos_pch <- c(rep(16, 4), rep(16, 2))
pt_cex <- 0.5
state_colors <- c(rgb(230, 159, 0, maxColorValue = 255), 
                  rgb(87, 180, 233, maxColorValue = 255))
CTDS_colors <- c("#a50026", "#313695")
## get CA border and save for future use
# library(tigris)
# options(tigris_use_cache = TRUE)
# US <- states(cb = T)
# CA <- US[US$NAME == "California", ]
# save(CA, file = "../data/whale/CA.RData")
load("../data/whale/CA.RData")
CA_proj <- st_transform(CA, st_CRS)
# map_whale <- get_googlemap(center = c(-121.1, 34.7), zoom = 7, maptype = "terrain", col = "bw")
# save(map_whale, file = "../data/whale/map.RData")
load("../data/whale/map.RData")
## load blue + fin whale data ----
whale <- read.csv("../data/whale/Blue and fin whales Southern California 2014-2015 - Argos data.csv")
whale <- whale[whale$argos.lc != "", ]
whale$POSIX <- as.POSIXct(whale$timestamp)
whale$argos.lc <- as.factor(whale$argos.lc)
ids <- unique(whale$individual.local.identifier)
## project data ----
coordinates(whale) <- c('location.long', 'location.lat')
proj4string(whale) <- CRS("+proj=longlat")
whale_st <- st_as_sf(whale)
whale_proj <- st_transform(whale_st, st_CRS)
## select individual ----
id <- ids[13]
time_ind <- 1:319 # 13 goes off DEM raster at the end of study period
whale_proj_id <- whale_proj[whale_proj$individual.local.identifier == id, ][time_ind, ]
## depth covariate ----
## data at: https://doi.org/10.7289/V500001J
# nc_data <- nc_open("../data/whale/crm_vol6.nc")
# dem <- raster("../data/whale/crm_vol6.nc")
# crs(dem) <- CRS("+proj=longlat")
# dem_proj <- projectRaster(dem, crs = st_CRS)
# dem_agg <- aggregate(dem_proj, 80) ## for CTDS model
# save(dem_proj, file = "../data/whale/dem_proj.RData")
# save(dem_agg, file = "../data/whale/dem_agg.RData")
# nc_close(nc_data)
load("../data/whale/dem_proj.RData")
load("../data/whale/dem_agg.RData") ## for CTDS model
## device ----
pdf("../fig/AMS_Notices/whale_obs.pdf")
## plot observed trajectory ----
## choose individual
whale_id <- whale[whale$individual.local.identifier == id, ][time_ind, ]
ggmap(map_whale, extent = "device") + 
  geom_point(aes(x = location.long, y = location.lat), 
            data = as.data.frame(whale_id), show.legend = F, 
            pch = argos_pch[whale_id$argos.lc],
            col = argos_colors[whale_id$argos.lc]) +
  geom_path(aes(x = location.long, y = location.lat), size = 0.5, lty = 3,
             data = as.data.frame(whale_id), show.legend = F, col = colors['obs']) +
  geom_point(aes(x = location.long, y = location.lat), 
             data = as.data.frame(whale_id)[c(1, nrow(whale_id)), ],
             col = colors['obs'], size = 4, shape = c(19, 17))
## dev.off ----
dev.off()
## functional ----
t <- as.numeric(whale_proj_id$POSIX)
pred_times <- data.frame(t = as.numeric(seq(min(whale_proj_id$POSIX), max(whale_proj_id$POSIX), "hour")))
s <- st_coordinates(whale_proj_id)
df <- data.frame(s_1 = s[, 1], s_2 = s[, 2], t = t)
fit_gam_r <- gam(list(s_1 ~ s(t, k = 85, fx = T, bs = "cr"), 
                      s_2 ~ s(t, k = 85, fx = T, bs = "cr")), 
                 family = mvn(d = 2), data = df)
r <- predict(fit_gam_r, pred_times, se.fit = T)
Sigma_s <- solve(crossprod(fit_gam_r$family$data$R))
## device ----
pdf("../fig/AMS_Notices/whale_gam.pdf")
## plot GAM fit ----
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(s, asp = 1, type = "n", xlab = "", ylab = "", axes = F, xlim = c(-60, 420))
rect(par()$usr[1], par()$usr[3], par()$usr[2], par()$usr[4], col = "gray90", border = F)
plot(st_geometry(CA_proj), add = T, col = "white")
points(s[c(1, nrow(s)), ], pch = c(17, 16), cex = 1.5)
lines(r$fit, col = colors['true'])
for(t in 1:nrow(pred_times)){
  lines(ellipse::ellipse(diag(r$se.fit[t, ]), centre = r$fit[t, ]),
        col = scales::alpha(colors['true'], 0.5))
}
points(s, cex = pt_cex, pch = argos_pch[whale_proj_id$argos.lc], 
       col = argos_colors[whale_proj_id$argos.lc])
legend("bottomleft", pch = c(16, 17, 16), pt.cex = c(pt_cex, 1.5, 1.5),
       legend = c("observed locations", "start", "stop"), bty = "n")
## compass + scale + cities
scale_bb <- c(100, 3640, 200, 3645)
rect(xleft = scale_bb[1], ybottom = scale_bb[2], xright = scale_bb[3], ytop = scale_bb[4])
rect(xleft = scale_bb[1] + 25, ybottom = scale_bb[2], xright = scale_bb[1] + 50, ytop = scale_bb[4], col = "black")
rect(xleft = scale_bb[1] + 75, ybottom = scale_bb[2], xright = scale_bb[3], ytop = scale_bb[4], col = "black")
text(c(scale_bb[1] + 50, scale_bb[3] + 7), rep(3654, 2), labels = c("50", "100km"))
arrows(x0 = -63, y0 = 4080, y1 = 4120, lwd = 2, length = 0.14)
text(-63, 4133, labels = "N", cex = 1.5)
points(c(388.98060, 62.83208), c(3768.420, 4145.279), pch = 15, cex = 1.2)
text(c(388.98060, 62.83208) + 12, c(3768.420, 4145.279) - 10, labels = c("Los Angeles", "San Jose"))
##
par(new = T, fig = c(0.45, 1, 0.8, 1))
plot(whale_proj_id$POSIX, st_coordinates(whale_proj_id)[, 1], pch = 16, cex = pt_cex, 
     col = argos_colors[whale_proj_id$argos.lc], xlab = "", 
     ylab = "", bty = "n", axes = F)
mtext("easting", 2, 2.2)
axis(2, at = c(0, 200) + par()$usr[3], labels = c(0, 200))
lines(pred_times$t, r$fit[, 1], col = colors['true'])
lines(pred_times$t, r$fit[, 1] - 1.96 * r$se.fit[, 1], lty = 2, col = scales::alpha(colors['true'], 0.5))
lines(pred_times$t, r$fit[, 1] + 1.96 * r$se.fit[, 1], lty = 2, col = scales::alpha(colors['true'], 0.5))
par(new = T, fig = c(0.45, 1, 0.6, 0.8))
plot(whale_proj_id$POSIX, st_coordinates(whale_proj_id)[, 2], pch = 16, cex = pt_cex, 
     col = argos_colors[whale_proj_id$argos.lc], xlab = "", yaxt = "n",
     ylab = "", bty = "n")
mtext("northing", 2, 2.2)
mtext("time", 1, 2.2)
axis(2, at = c(0, 200, 400) + par()$usr[3], labels = c(0, 200, 400))
lines(pred_times$t, r$fit[, 2], col = colors['true'])
lines(pred_times$t, r$fit[, 2] - 1.96 * r$se.fit[, 2], lty = 2, col = scales::alpha(colors['true'], 0.5))
lines(pred_times$t, r$fit[, 2] + 1.96 * r$se.fit[, 2], lty = 2, col = scales::alpha(colors['true'], 0.5))
## dev.off ----
dev.off()
## crawl ----
err_model = list(x = ~ 1)
err_model = list(x = ~ argos.lc - 1)
displayPar(err.model = err_model, data = whale_proj_id, Time.name = 'POSIX')
crwMLE_whales <- crwMLE(err.model = err_model, 
                        data = whale_proj_id, Time.name = 'POSIX', 
                        fixPar = rep(NA, 8), 
                        theta = log(c(2, 1.5, 0.5, 0.25, 5, 5, 5, 1)),
                        need.hess = T)
predTime <- seq(min(whale_proj_id$POSIX), max(whale_proj_id$POSIX), by = "1 hour")
n_time_pts <- length(predTime)
K <- 32
sim_whales_obj <- crwSimulator(crwMLE_whales, predTime = predTime)
sim_mu <- array(sapply(1:K, function(k){
  sim_whales <- crwPostIS(sim_whales_obj, fullPost = T)
  sim_whales$alpha.sim[c(1, which(sim_whales$locType == "p")), c('mu.x', 'mu.y')]
}), dim = c(n_time_pts, 2, K))
## table of CIs for argos variances ----
sapply(c(3:0, "A", "B"), function(class){
  mean(whale_proj_id$argos.error.radius[whale_proj_id$argos.lc == class])
}) * 1e-3
sqrt(exp(crwMLE_whales$par[c(4:1, 5:7)]))
sqrt(t(exp(crwMLE_whales$ci[c(4:1, 5:7), ])))
exp(crwMLE_whales$par[8])
exp(crwMLE_whales$ci[8, ])
## dev.off ----
dev.off()
## HMM priors + constants ----
n_hidden_states <- 2
stepPar0 <- c(2, 7, 2, 7)
anglePar0 <- rep(1, n_hidden_states)
## single-stage ----
k <- 2
elev_sim <- scale(extract(dem_proj, sim_mu[, , k]), center = F)
df <- data.frame(x = sim_mu[, 1, k], y = sim_mu[, 2, k], 
                 time = predTime, elev = elev_sim)
pred_df <- prepData(df, type = "UTM")
fit_k <- fitHMM(data = pred_df,                  
                nbStates = n_hidden_states,
                formula = ~elev,
                retryFits = 0,
                dist = list(step = "gamma", angle = "vm"), 
                Par0 = list(step = stepPar0, angle = anglePar0))
# plot(fit_k)
## HMM two-stage ----
cl <- makeForkCluster(8)
dem_proj_scaled <- scale(dem_proj, center = F)
MI_data <- parSapply(cl, 1:K, function(k){
  elev_sim <- extract(dem_proj, sim_mu[, , k]) / 1e3
  df <- data.frame(x = sim_mu[, 1, k], y = sim_mu[, 2, k], 
                   time = predTime, elev = elev_sim)
  prepData(df, type = "UTM")
}, simplify = F) 
stopCluster(cl)
## check to make sure paths didn't leave DEM extent
which(unlist(lapply(MI_data, function(x) sum(is.na(x$elev)))) > 0)
system.time({
  MIfit <- MIfitHMM(miData = MI_data, nSims = K, ncores = 8,
                    nbStates = n_hidden_states, poolEstimates = F,
                    formula = ~elev, retryFits = 1, 
                    dist = list(step = "gamma", angle = "vm"), 
                    Par0 = list(step = stepPar0, angle = anglePar0))
})
MIfit_pooled <- MIpool(MIfit, ncores = 4)
## device ----
pdf("../fig/AMS_Notices/HMM_summary.pdf", width = 5, height = 5)
## summary plots ----
plot(MIfit_pooled, plotCI = T, ask = F)
## dev.off ----
dev.off()
## device ----
pdf("../fig/AMS_Notices/whale_crawl.pdf")
## plot crawl fit w/HMM ----
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(s, asp = 1, type = "n", xlab = "", ylab = "", axes = F,
     xlim = c(-60, 420))
rect(par()$usr[1], par()$usr[3], par()$usr[2], par()$usr[4], col = "gray90", border = F)
plot(st_geometry(CA_proj), add = T, col = "white")
for(k in 1:K){
  lines(sim_mu[, , k], col = scales::alpha(colors['true'], 4/K))
}
points(s, cex = pt_cex, pch = argos_pch[whale_proj_id$argos.lc], 
       col = argos_colors[whale_proj_id$argos.lc])
points(s[c(1, nrow(s)), ], pch = c(17, 16), cex = 1.5)
for(t in 2:nrow(pred_times)){
  segments(x0 = MIfit_pooled$data[t - 1, 'x'], x1 = MIfit_pooled$data[t, 'x'], 
           y0 = MIfit_pooled$data[t - 1, 'y'], y1 = MIfit_pooled$data[t, 'y'], 
           col = state_colors[MIfit_pooled$Par$states[t - 1]], lwd = 2)
}
legend("bottomleft", pch = c(NA, NA, 16, 17, 16), pt.cex = c(NA, NA, pt_cex, 1.5, 1.5), 
       lty = c(1, 1, NA, NA, NA), lwd = c(1.5, 1.5, NA, NA, NA), col = c(state_colors, 1, 1, 1),
       legend = c("state 1", "state 2", "observed locations", "start", "stop"), bty = "n")
## compass + scale + cities
scale_bb <- c(100, 3640, 200, 3645)
rect(xleft = scale_bb[1], ybottom = scale_bb[2], xright = scale_bb[3], ytop = scale_bb[4])
rect(xleft = scale_bb[1] + 25, ybottom = scale_bb[2], xright = scale_bb[1] + 50, ytop = scale_bb[4], col = "black")
rect(xleft = scale_bb[1] + 75, ybottom = scale_bb[2], xright = scale_bb[3], ytop = scale_bb[4], col = "black")
text(c(scale_bb[1] + 50, scale_bb[3] + 7), rep(3654, 2), labels = c("50", "100km"))
arrows(x0 = -63, y0 = 4080, y1 = 4120, lwd = 2, length = 0.14)
text(-63, 4133, labels = "N", cex = 1.5)
points(c(388.98060, 62.83208), c(3768.420, 4145.279), pch = 15, cex = 1.2)
text(c(388.98060, 62.83208) + 12, c(3768.420, 4145.279) - 10, labels = c("Los Angeles", "San Jose"))
##
par(new = T, fig = c(0.45, 1, 0.8, 1))
plot(whale_proj_id$POSIX, st_coordinates(whale_proj_id)[, 1], pch = 16, cex = pt_cex, 
     col = argos_colors[whale_proj_id$argos.lc], xlab = "", 
     ylab = "", bty = "n", axes = F)
mtext("easting", 2, 2.2)
axis(2, at = c(0, 200) + par()$usr[3], labels = c(0, 200))
for(k in 1:K){
  lines(predTime, sim_mu[, 1, k], col = scales::alpha(colors['true'], 1/K))
}
par(new = T, fig = c(0.45, 1, 0.6, 0.8))
plot(whale_proj_id$POSIX, st_coordinates(whale_proj_id)[, 2], pch = 16, cex = pt_cex, 
     col = argos_colors[whale_proj_id$argos.lc], xlab = "", yaxt = "n",
     ylab = "", bty = "n")
mtext("northing", 2, 2.2)
mtext("time", 1, 2.2)
axis(2, at = c(0, 200, 400) + par()$usr[3], labels = c(0, 200, 400))
for(k in 1:K){
  lines(predTime, sim_mu[, 2, k], 
        col = scales::alpha(colors['true'], 1/K))
}
## dev.off ----
dev.off()
## ----
## CTDS ----
get_cells <- function(r, raster, t, n_out = 1e6){
  tout <- seq(min(t), max(t), l = n_out)
  r_interp <- apply(r, 2, function(r_c){
    approx(x = t, y = r_c, xout = tout)$y
  })
  extract(raster, r_interp, cellnumbers = T)
}
whale_DS1 <- get_cells(sim_mu[, , 1], dem_agg, predTime)
whale_DS2 <- get_cells(sim_mu[, , 2], dem_agg, predTime)
whale_cell_seq <- stack(dem_agg, dem_agg)
whale_cell_seq <- setValues(whale_cell_seq, NA)
whale_cell_seq[[1]][whale_DS1[, 1]] <- 1
whale_cell_seq[[2]][whale_DS2[, 1]] <- 1
## device ----
png("../fig/AMS_Notices/whale_ctds.png", width = 960, height = 960)
## plot example ---
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(s, asp = 1, type = "n", xlab = "", ylab = "", axes = F, 
     xlim = c(115, 227), ylim = c(3670, 3900))
plot(dem_agg, add = T, legend = F, col = gray.colors(1e2, 0, 1, alpha = 0.75))
plot(st_geometry(CA_proj), add = T, col = "white")
plot(whale_cell_seq[[1]], col = scales::alpha(CTDS_colors[1], 0.5), add = T, legend = F)
plot(whale_cell_seq[[2]], col = scales::alpha(CTDS_colors[2], 0.5), add = T, legend = F)
lines(sim_mu[, , 1], col = CTDS_colors[1], lwd = 3)
lines(sim_mu[, , 2], col = CTDS_colors[2], lwd = 3)
## dev.off ----
dev.off()
