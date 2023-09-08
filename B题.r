library(measurements)

# 第一问
depth_delta <- tan(pi / 120) * 200
raw_result1 <- read.csv("result1.csv", fileEncoding = "UTF-8-BOM", row.names = 1)
raw_result1$X.800 <- c(NA, NA, 0)
raw_result1[1, ] <- 4:-4 * depth_delta + 70
raw_result1[2, ] <- raw_result1[1, ] * (sinpi(5 / 6) * tan(pi / 3) / sinpi(19 / 120) + sin(pi / 3) / sinpi(7 / 40))
raw_result1[3, ] <- 1 - 200 / raw_result1[2, ]
write.csv(raw_result1, "result1_raw.csv")
raw_result1[2, ] <- raw_result1[2, ] * cos(pi / 120)
raw_result1[3, ] <- 1 - 200 / raw_result1[2, ]
write.csv(raw_result1, "result1_cos.csv")
# 第一问画图
plot_delta <- tan(pi / 3) * 100
plot(c(-800, 800), c(-90, 0), "n", cex.axis = 2, cex.lab = 1.1, las = 1, xlab = "海底", ylab = "海拔")
for (i in seq(-800, 800, 200)) {
    polygon(c(i - plot_delta, i, i + plot_delta), c(-100, 0, -100), col = rgb(0, 0, 0, 0.1), lty = "dotted")
}
lines(c(-1000, 1000), c(0, 0), lty = "dashed")
polygon(c(-1000, 1000, 1000, -1000), c(-(depth_delta * 5 + 70), depth_delta * 5 - 70, -100, -100), col = rgb(1, 1, 1, 0.75))


# 第二问
raw_result2 <- read.csv("result2.csv", fileEncoding = "UTF-8-BOM")
colnames(raw_result2) <- raw_result2[1, ]
named_result2 <- raw_result2[-1, ]
row.names(named_result2) <- named_result2[, 2]
matrix_result2 <- named_result2[, -(1:2)]
depth_cos <- cospi(as.integer(row.names(matrix_result2)) / 180)
depth_angle <- atan(depth_cos * tan(pi / 120))
cover_result2 <- (depth_cos %*% t(conv_unit(as.double(colnames(matrix_result2)), "naut_mi", "m")) * tan(pi / 120) + 120) * (sinpi(5 / 6) * tan(pi / 3) / sin(pi / 6 - depth_angle) + sin(pi / 3) / sinpi(pi / 6 + depth_angle))
write.csv(cover_result2, "result2_raw.csv")


raw_attachment <- read.csv("附件.csv", fileEncoding = "UTF-8-BOM")
colnames(raw_attachment) <- raw_attachment[1, ]
raw_attachment <- raw_attachment[-1, ]
row.names(raw_attachment) <- raw_attachment[, 2]
raw_attachment <- raw_attachment[, -(1:2)]
contour(pin = c(1, 1), z = as.matrix(raw_attachment))
