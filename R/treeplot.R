plot.tree = function(phy, ...) {
    vargs = list(...)
    method = ifelse(is.null(method <- vargs$method), "cartesian",
        match.arg(method, c("cartesian", "polar")))
    if (is.null(edge.color <- vargs$edge.color)) edge.color = rep(1, Nnode(phy))
    if (is.null(edge.width <- vargs$lwd)) edge.width = rep(1, Nnode(phy))
    if (length(edge.color) < Nnode(phy)) edge.color = rep(edge.color[1], Nnode(phy))
    if (length(edge.width) < Nnode(phy)) edge.width = rep(edge.width[1], Nnode(phy))
    if (is.null(tip.labels <- vargs$tip.labels)) tip.labels = FALSE
    if (is.null(cex.label <- vargs$cex.label)) cex.label = 1
    if (method == "cartesian") {
        if (is.null(direction <- vargs$direction)) direction = "right"
        direction = match.arg(direction, c("up", "down", "left", "right"))
        xlim = switch(direction,
            up = c(0, Ntip(phy)+1),
            down = c(0, Ntip(phy)+1),
            left = c(min(ages(phy)), max(ages(phy))),
            right = c(min(ages(phy)), max(ages(phy))))
        ylim = switch(direction,
            up = c(min(ages(phy)), max(ages(phy))),
            down = c(min(ages(phy)), max(ages(phy))),
            left = c(0, Ntip(phy)+1),
            right = c(0, Ntip(phy)+1))
        dir = switch(direction,
            up = 2L,
            down = 3L,
            left = 1L,
            right = 0L)
        L = .Call(rcl_plot_ctree, phy, ages(phy), dir)
        edge = L[[1L]]
        bar = L[[2L]]
        plot.new()
        plot.window(xlim=xlim, ylim=ylim)
        segments(edge[, 1L], edge[, 3L], edge[, 2L], edge[, 4L],
            lwd=edge.width, col=edge.color)
        segments(bar[, 1L], bar[, 3L], bar[, 2L], bar[, 4L],
            lwd=edge.width[-(1:Ntip(phy))], col=edge.color[-(1:Ntip(phy))])
        if (tip.labels) {
            pos = switch(direction,
                up = 3L,
                down = 1L,
                left = 2L,
                right = 4L)
            text(edge[1L:Ntip(phy), 1L], edge[1L:Ntip(phy), 3L], labels=tiplabels(phy), cex=cex.label, pos=pos)
        }
        invisible(L)
    } else {
        vtheta = ifelse(is.null(vtheta <- vargs$vtheta), 5 * (pi/180), vtheta * (pi/180))
        rbf = ifelse(is.null(rbf <- vargs$rbf), 0.01, rbf)
        vstep = (2*pi - vtheta) / (Ntip(phy) - 1)
        theta = .Call(rcl_plot_ptree, phy, vstep)
        tree.height = max(ages(phy))
        rb = rbf * tree.height
        plot.new()
        plot.window(xlim=c(-(tree.height + rb), tree.height+rb), ylim=c(-(tree.height+rb), tree.height+rb), asp=1)
        x0 = (rb + ages(phy)) * cos(theta[, 1L])
        y0 = (rb + ages(phy)) * sin(theta[, 1L])
        x1 = (ages(phy) + rb - brlens(phy)) * cos(theta[, 1L])
        y1 = (ages(phy) + rb - brlens(phy)) * sin(theta[, 1L])
        segments(x0, y0, x1, y1, lwd=edge.width, col=edge.color)
        draw.arc(ages(phy)[-(1:Ntip(phy))] + rb, theta[-(1:Ntip(phy)), 2:3],
            lwd=edge.width[-(1:Ntip(phy))], col=edge.color[-(1:Ntip(phy))])
        invisible(theta)
    }
}


draw.arc = function(r, theta, grain=30, ...) {
    path = t(apply(theta, 1, function(z) seq(z[1], z[2], length.out=grain+1)))
    arcs = apply(cbind(r, theta), 1, function(z) {
        zz = seq(z[2], z[3], length.out=grain+1)
        arcs = cbind(zz[-(grain+1)], zz[-1])
        x0 = z[1] * cos(arcs[, 1])
        y0 = z[1] * sin(arcs[, 1])
        x1 = z[1] * cos(arcs[, 2])
        y1 = z[1] * sin(arcs[, 2])
        segments(x0, y0, x1, y1, ...)
    })
}
