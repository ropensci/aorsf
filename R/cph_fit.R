

cph_fit <- function(x,
                    y,
                    weights,
                    x_transforms,
                    method = 1,
                    eps = 1e-09,
                    iter_max = 20,
                    rescale = TRUE){

        x_transforms = x_scale_wtd(x, weights)

        newtraph_cph(x,
                     y,
                     weights,
                     x_transforms,
                     method,
                     eps,
                     iter_max,
                     rescale)

}
