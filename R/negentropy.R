negentropy <-
function(x){(1/12)*(mean(x^3)^2)+(1/48)*kurtosis(x)^2}
