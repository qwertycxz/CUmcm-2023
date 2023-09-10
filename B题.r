library(measurements)
library(nloptr)

# 第一问
depth_delta <- tan(pi / 120) * 200
raw_result1 <- read.csv("result1.csv", fileEncoding = "UTF-8-BOM", row.names = 1)
raw_result1$X.800 <- c(NA, NA, 0)
raw_result1[1, ] <- 4:-4 * depth_delta + 70
raw_result1[2, ] <- raw_result1[1, ] * (sinpi(5 / 6) * tan(pi / 3) / sinpi(19 / 120) + sin(pi / 3) / sinpi(7 / 40))
raw_result1[3, 2:9] <- 1 - sinpi(5 / 6) * 200 / raw_result1[2, ] / sinpi(19 / 120)
write.csv(raw_result1, "result1_raw.csv")
raw_result1[2, ] <- raw_result1[2, ] * cos(pi / 120)
raw_result1[3, 2:9] <- 1 - sinpi(5 / 6) * 200 / raw_result1[2, ] / sinpi(19 / 120)
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
# 第二问画图
plot_delta <- conv_unit(tan(pi / 3) * 250, "m", "naut_mi")
plot_left <- -depth_cos * conv_unit(-0.3, "naut_mi", "m") * tan(pi / 120) - 120
plot_right <- -depth_cos * conv_unit(2.4, "naut_mi", "m") * tan(pi / 120) - 120
plot_color <- rgb(t(col2rgb(c("red", "orange", "yellow", "green", "blue")) / 255), alpha = 0.25)
plot(c(0, 2.1), c(-220, 0), "n", cex.axis = 2, cex.lab = 1.1, las = 1, xlab = "海底", ylab = "海拔")
for (i in 5:1) {
    polygon(c(-0.3, 2.4, 2.4, -0.3), c(plot_left[i], plot_right[i], -250, -250), col = plot_color[i])
}
for (i in seq(0, 2.1, 0.3)) {
    polygon(c(i - plot_delta, i, i + plot_delta), c(-250, 0, -250), col = rgb(0, 0, 0, 0.1), lty = "dotted")
}
lines(c(-1000, 1000), c(0, 0), lty = "dashed")


# 第三问
# 最简情况
size_x <- conv_unit(4, "naut_mi", "m")
size_y <- conv_unit(2, "naut_mi", "m")
depth_delta <- tan(pi / 120)
cover_lines <- NULL
pos_x <- tan(pi / 3) * (depth_delta * size_x / 2 + 110) - size_x / 2
x_left <- pos_x
while (pos_x * 2 < size_x) {
    cover_depth <- -depth_delta * pos_x + 110
    cover_left <- cover_depth * sinpi(5 / 6) * tan(pi / 3) / sinpi(19 / 120)
    cover_right <- cover_depth * sin(pi / 3) / sinpi(7 / 40)
    cover_distance <- (cover_left + cover_right) * sinpi(19 / 120) * 0.9 / sinpi(5 / 6)
    cover_lines <- c(cover_lines, pos_x)
    pos_x <- pos_x + cover_distance
}
if ((cover_right + cover_lines[length(cover_lines)]) * 2 < size_x) {
    cover_lines <- c(cover_lines, size_x / 2)
}
# 画图
plot_delta_x <- conv_unit(tan(pi / 3) * 250, "m", "naut_mi")
plot_delta_y <- conv_unit(depth_delta * 4, "naut_mi", "m")
plot(c(-2, 2), c(-200, 0), "n", cex.axis = 2, cex.lab = 1.1, las = 1, xlab = "海底", ylab = "")
for (i in conv_unit(cover_lines, "m", "naut_mi")) {
    polygon(c(i - plot_delta_x, i, i + plot_delta_x), c(-250, 0, -250), col = rgb(0, 0, 0, 0.1), lty = "dotted")
}
polygon(c(-4, 4, 4, -4), c(-(plot_delta_y + 110), plot_delta_y - 110, -250, -250), col = rgb(1, 1, 1, 0.75))
lines(c(-2, -2, 2, 2), c(-250, 0, 0, -250), lty = "dashed")
# 一般情况
size_x <- conv_unit(4, "naut_mi", "m")
size_y <- conv_unit(2, "naut_mi", "m")
getDensestLine <- function(angle, as_goal = TRUE) {
    # size_x_transformed <- size_x / sin(angle)
    size_x_prolong <- (size_x + size_y / tan(angle)) / 2
    depth_delta <- cos(angle - pi / 2) * tan(pi / 120)
    depth_sin1 <- sin(pi / 6 - atan(depth_delta))
    depth_sin2 <- sin(pi / 6 + atan(depth_delta))
    cover_lines <- NULL
    cover_distances <- NULL
    pos_x <- (size_x / 2 + (tan(pi / 3) * (depth_delta * size_x / 2 + 110) - size_x / 2)) / sin(angle) - size_x / 2
    while (pos_x < size_x_prolong) {
        cover_depth <- -depth_delta * pos_x * sin(angle) + 110
        cover_left <- cover_depth * sinpi(5 / 6) * tan(pi / 3) / depth_sin1
        cover_right <- cover_depth * sin(pi / 3) / depth_sin2
        cover_distance <- (cover_left + cover_right) * depth_sin1 * 0.9 / sinpi(5 / 6)
        cover_distances <- c(cover_distances, cover_distance)
        cover_lines <- c(cover_lines, pos_x)
        pos_x <- pos_x + cover_distance * cos(angle - pi / 2) / sin(angle)
    }
    if (cover_right + cover_lines[length(cover_lines)] < sin(angle) * size_x_prolong) {
        cover_lines <- c(cover_lines, size_x_prolong)
        cover_distances <- c(cover_distances, size_x / sin(angle) - cover_lines[length(cover_lines)])
    }
    if (as_goal) {
        return(cover_lines)
    }
    return(cover_distances)
}
cover_lines <- getDensestLine(pi / 2)
cover_distances <- getDensestLine(pi / 2, FALSE)
write.csv(cover_distances, "result3_raw.csv")
getDensestPosition <- function(angle) {
    cover_lines <- getDensestLine(angle)
    position <- matrix(1:4, 1)[-1, ]
    colnames(position) <- c("x1", "y1", "x2", "y2")
    for (i in cover_lines) {
        x_pos <- NULL
        y_pos <- NULL
        x_bottom <- i - size_y / tan(angle)
        y_left <- size_y / 2 - (size_x / 2 + i) * tan(angle)
        y_right <-  (size_x / 2 - i) * tan(angle) + size_y / 2
        if (i > -size_x / 2 && i < size_x / 2) {
            x_pos <- c(x_pos, i)
            y_pos <- c(y_pos, size_y / 2)
        }
        if (x_bottom > -size_x / 2 && x_bottom < size_x / 2) {
            x_pos <- c(x_pos, x_bottom)
            y_pos <- c(y_pos, -size_y / 2)
        }
        if (y_right > -size_y / 2 && y_right < size_y / 2) {
            x_pos <- c(x_pos, size_x / 2)
            y_pos <- c(y_pos, y_right)
        }
        if (y_left > -size_y / 2 && y_left < size_y / 2) {
            x_pos <- c(x_pos, -size_x / 2)
            y_pos <- c(y_pos, y_left)
        }
        position <- rbind(position, c(x_pos[1], y_pos[1], x_pos[2], y_pos[2]))
    }
    return(position)
}
getDensestPosition(pi / 2)
getTotalLength <- function(angle) {
    cover_position <- getDensestPosition(angle)
    return(sum(sqrt(apply(cover_position[, c(1, 3)], 1, diff) ^ 2 + apply(cover_position[, c(2, 4)], 1, diff) ^ 2)))
}
getTotalLength(pi / 2)
getCoincidenceRate <- function(angle) {
    depth_delta <- cos(angle - pi / 2) * tan(pi / 120)
    depth_sin1 <- sin(pi / 6 - atan(depth_delta))
    depth_sin2 <- sin(pi / 6 + atan(depth_delta))
    cover_position <- getDensestPosition(angle)
    depth_delta <- tan(pi / 120)
    x_deepest <- apply(cover_position[, c(1, 3)], 1, min)
    depth_deepest <- -depth_delta * x_deepest + 110
    # x_swallowest <- apply(cover_position[, c(1, 3)], 1, max)
    # depth_swallowest <- -depth_delta * x_swallowest + 110
    # return(depth_deepest * 0.1 / depth_swallowest)
    cover_left <- depth_deepest * sinpi(5 / 6) * tan(pi / 3) / depth_sin1
    cover_right <- depth_deepest * sin(pi / 3) / depth_sin2
    cover_distances <- getDensestLine(angle, FALSE)
    return(1 - sinpi(5 / 6) * cover_distances / (cover_left + cover_right) / depth_sin1)
}
getCoincidenceRate(pi / 2)
# NLopt
evalGIneq <- function(angle) {
    return(max(getCoincidenceRate(angle)) - 0.2)
}
nlopt <- nloptr(pi / 2, getTotalLength, eval_g_ineq = evalGIneq, lb = 0, ub = pi / 2, opts = list("algorithm" = "NLOPT_GN_ORIG_DIRECT"))
nlopt <- nloptr(pi / 2, getTotalLength, eval_g_ineq = evalGIneq, lb = 0, ub = pi / 2, opts = list("algorithm" = "NLOPT_GN_ESCH"))
nlopt
length_data <- matrix(1:2, 1)[-1, ]
for (i in seq(pi / 180, pi / 2, length.out = 1000)) {
    length_data <- rbind(length_data, c(i * 180 / pi, getTotalLength(i)))
}
plot(length_data, cex.axis = 2, cex.lab = 1.1, type = "l", xlab = "β", ylab = "总长度")
length_data <- matrix(1:2, 1)[-1, ]
for (i in seq(pi / 4, pi / 2, length.out = 1000)) {
    length_data <- rbind(length_data, c(i * 180 / pi, getTotalLength(i)))
}
plot(length_data, cex.axis = 2, cex.lab = 1.1, xlab = "β", ylab = "总长度")
# fun <- function(angle) {
#     depth_cos <- cos(angle)
#     depth_angle <- atan(depth_cos * tan(pi / 120))
#     (depth_cos * conv_unit(2.1, "naut_mi", "m") * tan(pi / 120) + 120) * (sinpi(5 / 6) * tan(pi / 3) / sin(pi / 6 - depth_angle) + sin(pi / 3) / sinpi(pi / 6 + depth_angle))
# }


# 第四问
raw_attachment <- read.csv("附件.csv", fileEncoding = "UTF-8-BOM")
colnames(raw_attachment) <- raw_attachment[1, ]
named_attachment <- raw_attachment[-1, ]
row.names(named_attachment) <- named_attachment[, 2]
matrix_attachment <- named_attachment[, -(1:2)]
contour(z = as.matrix(matrix_attachment))
getYDepth <- function(y) {
    pos <- NA
    if (y %in% row.names(matrix_attachment)) {
        pos <- matrix_attachment[as.character(y), ]
    }
    else {
        y_floor <- floor(y * 50) / 50
        pos <- (matrix_attachment[as.character(y_floor + 0.02), ] - matrix_attachment[as.character(y_floor), ]) * (y - y_floor) / 0.02 + matrix_attachment[as.character(y_floor), ]
    }
    return(pos)
}
getDepth <- function(x, y = NULL) {
    if (is.null(y)) {
        y <- x[2]
        x <- x[1]
    }
    pos <- getYDepth(y)
    if (x %in% colnames(matrix_attachment)) {
        return(pos[as.character(x)])
    }
    x_floor <- floor(x * 50) / 50
    return((pos[as.character(x_floor + 0.02)] - pos[as.character(x_floor)]) * (x - x_floor) / 0.02 + pos[as.character(x_floor)])
}
# 东西向航线(北深南浅)
size_x <- conv_unit(4, "naut_mi", "m")
size_y <- conv_unit(5, "naut_mi", "m")
cover_lines <- NULL
cover_length <- 0
pos_y <- 5 # 海里
while (pos_y >= 0) {
    pos_x <- ceiling(pos_y * 40) # 序号(从0开始)
    cover_depth <- mean(unlist(getYDepth(pos_y)[1:pos_x + 1])) # 米
    depth_delta <- abs(cover_depth - 24.4) / size_y # 米
    depth_angle <- atan(depth_delta) # 弧度
    cover_left <- cover_depth * sinpi(5 / 6) * tan(pi / 3) / sin(pi / 6 - depth_angle) # 米
    cover_right <- cover_depth * sin(pi / 3) / sinpi(pi / 6 + depth_angle) # 米
    cover_distance <- (cover_left + cover_right) * sinpi(pi / 6 - depth_angle) * 0.9 / sinpi(5 / 6) # 米
    cover_lines <- c(cover_lines, pos_y) # 海里
    cover_length <- cover_length + pos_y * 0.8 # 海里
    pos_y <- pos_y - conv_unit(cover_distance, "m", "naut_mi")
}
cover_lines1 <- conv_unit(cover_lines, "naut_mi", "m")
cover_length1 <- conv_unit(cover_length, "naut_mi", "m")
if ((cover_right < cover_lines[length(cover_lines)])) {
    cover_lines1 <- c(cover_lines1, 0)
}
write.csv(cover_lines1, "result4_left_y_pos.csv")
# 西南-东北向航线(东南深西北浅)
cover_lines <- NULL
cover_length <- 0
pos_x <- 4 # 海里
while (pos_x >= 0) {
    seq_x <- seq(4, pos_x, -0.02) # 序列, 海里
    seq_y <- seq_x * 1.25 # 序列, 海里
    cover_depth <- mean(unlist(apply(cbind(seq_x, seq_y), 1, getDepth))) # 米
    depth_delta <- abs(cover_depth - 197.2) / sqrt((size_x / 2) ^ 2 + (size_y / 2) ^ 2) # 米
    depth_angle <- atan(depth_delta) # 弧度
    cover_left <- cover_depth * sinpi(5 / 6) * tan(pi / 3) / sin(pi / 6 - depth_angle) # 米
    cover_right <- cover_depth * sin(pi / 3) / sinpi(pi / 6 + depth_angle) # 米
    cover_distance <- (cover_left + cover_right) * sinpi(pi / 6 - depth_angle) * 0.9 / sinpi(5 / 6) # 米
    cover_lines <- c(cover_lines, pos_x) # 海里
    cover_length <- cover_length + pos_x * 0.8 # 海里
    pos_x <- pos_x - conv_unit(cover_distance, "m", "naut_mi")
}
cover_lines2 <- conv_unit(cover_lines, "naut_mi", "m")
cover_length2 <- conv_unit(cover_length, "naut_mi", "m")
if ((cover_right < cover_lines[length(cover_lines)])) {
    cover_lines2 <- c(cover_lines2, 0)
}
write.csv(cover_lines2[-length(cover_lines2)], "result4_right_x_pos.csv")
# 画图
plot_x <- NULL
add_zero <- TRUE
for (i in cover_lines1) {
    plot_x <- c(plot_x, i * 0.8)
    if (add_zero) {
        plot_x <- c(plot_x, 0, 0)
    }
    add_zero <- !add_zero
}
plot_xy_x <- NULL
add_zero <- TRUE
for (i in rev(cover_lines2)) {
    plot_xy_x <- c(plot_xy_x, i)
    if (add_zero) {
        plot_xy_x <- c(plot_xy_x, size_x, size_x)
    }
    add_zero <- !add_zero
}
plot_xy_y <- 0
add_zero <- FALSE
for (i in rev(cover_lines2)) {
    plot_xy_y <- c(plot_xy_y, i * 1.25)
    if (add_zero) {
        plot_xy_y <- c(plot_xy_y, 0, 0)
    }
    add_zero <- !add_zero
}
plot_xy_y <- rev(plot_xy_y)[-1]
plot(plot_x, rep(cover_lines1, each = 2), "l", cex.axis = 2, cex.lab = 1.1, xlab = "东西向(m)", ylab = "南北向(m)")
lines(c(size_x, plot_xy_x[-(1:2)]), c(size_y, plot_xy_y[-(1:2)]))
