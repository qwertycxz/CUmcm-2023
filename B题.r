# 第一问
depth_delta <- tan(pi / 120) * 200
raw_result1 <- read.csv("result1.csv", fileEncoding = "UTF-8-BOM", row.names = 1)
raw_result1$X.800 <- c(NA, NA, 0)
raw_result1[1, ] <- 4:-4 * depth_delta + 70
raw_result1[2, ] <- raw_result1[1, ] * cos(pi / 120) * sin(pi * 5 / 6) * tan(pi / 3) * 2 / sin(pi * 19 / 120)
raw_result1[3, ] <- 1 - 200 / raw_result1[2, ]
write.csv(raw_result1, "raw_result1.csv")

raw_result2 <- read.csv("result2.csv", fileEncoding = "UTF-8-BOM")
colnames(raw_result2) <- raw_result2[1, ]
raw_result2 <- raw_result2[-1, ]
row.names(raw_result2) <- raw_result2[, 2]
raw_result2 <- raw_result2[, -(1:2)]
raw_attachment <- read.csv("附件.csv", fileEncoding = "UTF-8-BOM")
colnames(raw_attachment) <- raw_attachment[1, ]
raw_attachment <- raw_attachment[-1, ]
row.names(raw_attachment) <- raw_attachment[, 2]
raw_attachment <- raw_attachment[, -(1:2)]
contour(pin = c(1, 1), z = as.matrix(raw_attachment))
