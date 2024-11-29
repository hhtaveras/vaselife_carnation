lmv.d.plot <- function(mapthis,
                       outfile,
                       mapthese = NULL,
                       at.axis = NULL,
                       autoconnadj = TRUE,
                       cex.axis = par("cex.axis"),
                       cex.lgtitle = par("cex.main"),
                       cex.main = par("cex.main"),
                       col.axis = par("col.axis"),
                       col.lgtitle = par("col.main"),
                       col.main = par("col.main"),
                       conndf = NULL,
                       denmap = FALSE,
                       dupnbr = FALSE,
                       font.axis = par("font.axis"),
                       font.lgtitle = par("font.main"),
                       font.main = par("font.main"),
                       header = TRUE,
                       labdist = .3,
                       labels.axis = TRUE,
                       lcex = par("cex"),
                       lcol = par("col"),
                       lfont = par("font"),
                       lgperrow = NULL,
                       lgtitles = NULL,
                       lgw = 0.25,
                       lg.col = NULL,
                       lg.lwd = par("lwd"),
                       lty.axis = "solid",
                       lwd.axis = 1,
                       lwd.ticks.axis = lwd.axis,
                       main = NULL,
                       markerformatlist = NULL,
                       maxnbrcolsfordups = 3,
                       pdf.bg = "transparent",
                       pdf.family = "Helvetica",
                       pdf.fg = "black",
                       pdf.width = NULL,
                       pdf.height = NULL,
                       pdf.pointsize = 12,
                       pdf.title = "LinkageMapView R output",
                       posonleft = NULL,
                       prtlgtitles = TRUE,
                       qtldf = NULL,
                       revthese = NULL,
                       rcex = par("cex"),
                       rcol = par("col"),
                       rfont = par("font"),
                       roundpos = 1,
                       rsegcol = TRUE,
                       ruler = FALSE,
                       sectcoldf = NULL,
                       segcol = NULL,
                       qtlscanone = NULL,
                       showonly = NULL,
                       units = "cM",
                       ylab = units
)

{
  pgx <- 0.5   # where on page x-axis to draw
  pgy <- 0.5   # where on page y-axis to draw
  
  oldpar <-
    par(no.readonly = TRUE)   # save parameters for restoring
  on.exit(par(oldpar))
  
  # ----------------- Begin edit arguments -----------------------------------
  if (!is.null(roundpos)) {
    if (!is.numeric(roundpos)) {
      stop("roundpos - number of places after decimal point - must be numeric")
    }
  }
  
  # edit and convert font from text if necessary
  
  font.axis <- convertfont("font.axis", font.axis)
  font.lgtitle <- convertfont("font.lgtitle", font.lgtitle)
  font.main <- convertfont("font.main", font.main)
  lfont <- convertfont("lfont", lfont)
  rfont <- convertfont("rfont", rfont)
  if (!is.null(markerformatlist$font)) {
    markerformatlist$font <-
      convertfont("markerformatlist$font", markerformatlist$font)
  }
  
  # read input for further edits ---
  
  if ("cross" %in% class(mapthis)) {
    if (is.null(mapthese)) {
      mapthese <- qtl::chrnames(mapthis)
    }
    lgin <- readlgcross(mapthis, mapthese)
  }
  else if ("character" %in% class(mapthis)) {
    if (is.null(mapthese)) {
      lgin <- readlgtext(mapthis, header = header)
      mapthese <- unique(lgin$group)
    } else{
      lgin <- readlgtext(mapthis, mapthese, header = header)
    }
  }
  else if ("data.frame" %in% class(mapthis)) {
    if (is.null(mapthese)) {
      lgin <- readlgdf(as.data.frame(mapthis))
      mapthese <- unique(lgin$group)
    } else{
      lgin <- readlgdf(mapthis, mapthese)
    }
  }
  else {
    stop("first parameter, mapthis, must be a data frame, a filename, or an r/qtl cross object")
  }
  
  if (!is.null(revthese)) {
    if (!all(revthese %in% mapthese)) {
      stop ("Requested chrnames to reverse not in those requested to plot")
    }
  }
  
  if (!is.null(lgperrow)) {
    if (!(is.numeric(lgperrow) && floor(lgperrow) == lgperrow)) {
      stop ("lgperrow must be NULL or an integer")
    }
  }
  else {
    lgperrow <- length(mapthese)
  }
  
  nbrrows <- ceiling(length(mapthese) / lgperrow)
  
  if (denmap) {
    autoconnadj <- FALSE
    rsegcol <- FALSE
    ruler <- TRUE
    showonly <- NULL
    markerformatlist <- NULL
    conndf <- NULL
  }
  
  # user can specify colors for segments
  if (!is.null(segcol)) {
    # make sure segcol column exists in input
    if (!(segcol %in% colnames(lgin))) {
      stop (c("segcol column ", segcol, " not found in mapthis"))
    }
    else {
      rsegcol <- FALSE
    }
  }
  
  if (!is.null(showonly)) {
    if (!all(showonly %in% lgin$locus)) {
      stop (c("Requested showonly locus not found:", showonly[!(showonly %in% lgin$locus)]))
    }
    else {
      # dups really screw up the code, correct for the user
      showonly <- unique(showonly)
    }
  }
  
  # else if input has null or NA Locus names convert existing locus names
  # to showonly list so position labels won't show for the null/NA
  # ones
  
  else {
    if (any(lgin$locus == "") || any(is.na(lgin$locus))) {
      notnull <- lgin$locus[which(lgin$locus != "")]
      showonly <- notnull[which(!is.na(notnull))]
    }
  }
  
  # make sure qtl data frame passed has no factors
  if (!is.null(qtldf)) {
    fas <- sapply(qtldf, is.factor)
    qtldf[fas] <- lapply(qtldf[fas], as.character)
  }
  # make sure qtlscanone is df and add it to (or create) qtldf
  if (!is.null(qtlscanone)) {
    if (!("data.frame" %in% class(qtlscanone) &
          "scanone" %in% class(qtlscanone))) {
      stop (c("qtlscanone should be a data.frame for r/qtl"))
    }
    else {
      qtldf <-
        usescanone(qtlscanone, qtldf, mapthese, pdf.fg, maxdec = roundpos)
    }
  }
  
  
  if (is.null(sectcoldf) && denmap == TRUE) {
    # use default density map coloring
    sectcoldf <- lmvdencolor(lgin)
  }
  # make sure sectcoldf data frame passed has no factors
  if (!is.null(sectcoldf)) {
    fas <- sapply(sectcoldf, is.factor)
    sectcoldf[fas] <- lapply(sectcoldf[fas], as.character)
  }
  
  # make sure vector of positioning requests is valid
  if (!is.null(posonleft)) {
    if (!(all(posonleft %in% c(T, F)) &&
          length(posonleft) == length(mapthese))) {
      stop(
        "Position on vector must be same length as # of linkage groups to map and only TRUE or FALSE"
      )
    }
  }
  else {
    posonleft <- rep(TRUE, length.out = length(mapthese))
  }
  
  # make sure ruler is valid
  if (!(is.null(ruler)) && !(length(ruler) == 1)) {
    stop("ruler must be a single TRUE or FALSE")
  }
  
  
  if (!(ruler == T) && !(ruler == F))              {
    stop("ruler may only be TRUE or FALSE")
  }
  
  # make sure dupnbr is valid
  if (!(is.null(dupnbr)) && !(length(dupnbr) == 1)) {
    stop("dupnbr must be a single TRUE or FALSE")
  }
  
  
  if (!(dupnbr == T) && !(dupnbr == F))              {
    stop("dupnbr may only be TRUE or FALSE")
  }
  
  if (dupnbr == T) {
    maxnbrcolsfordups <- 2
  }
  
  if (!(prtlgtitles == T) && !(prtlgtitles == F))              {
    stop("prtlgtitles may only be TRUE or FALSE")
  }
  # ----------------- End edit arguments -----------------------------------
  # ------------- Begin set par options --------------------
  
  # use top par("mar" for linkage group titles)
  # use top par("oma" for main title)
  # use left par("mar" for ruler)
  # use left par("oma" for ruler units)
  
  if (ruler) {
    leftmar <- 1 + ceiling(cex.axis)
    if (!is.null(units)) {
      leftmar <- leftmar + 1
    }
  }
  else {
    leftmar <- 1
  }
  
  if (prtlgtitles) {
    topmar <- ceiling(cex.lgtitle) + 1
  }
  else {
    topmar <- 1
  }
  
  if (!is.null(main)) {
    topoma <- 1+ ceiling(cex.main)
  }
  else {
    topoma <- 1
  }
  
  # ------------- End set par options --------------------
  
  pdf.options(
    bg = pdf.bg,
    title = pdf.title,
    family = pdf.family,
    pointsize = pdf.pointsize,
    fg = pdf.fg
  )
  on.exit(pdf.options(reset = TRUE), add = TRUE)
  
  # pdf size doesn't matter here - just for reqdim plotting
  # pdf will be reallocated before actually drawing
  
  pdf(outfile, width = 30, height = 30)
  on.exit(dev.off(), add = TRUE)
  
  par(mar = c(0, leftmar, topmar, 0),oma = c(1, 0, topoma, 1),new=FALSE)
  
  
  lg <- list()
  dim <- list()
  lrcol <- list()
  llcol <- list()
  lrfont <- list()
  llfont <- list()
  lrcex <- list()
  llcex <- list()
  qtl <- list()
  solist <- list()
  sectcol <- list()
  
  # make a list of linkage groups to process
  # reverse position numbers if requested
  for (i in 1:length(mapthese)) {
    lg[[i]] <- getlg(lgin, mapthese[i], dupnbr,roundpos)
    editlgdf(lg[[i]])  # display message and stop if not in correct format
    if (lg[[i]]$group[1] %in% revthese) {
      lg[[i]]$locus <- rev(lg[[i]]$locus)
      lg[[i]]$position <- revpos(lg[[i]]$position, roundpos)
      #reverse segcol column if it exists
      if (!is.null(segcol)) {
        lg[[i]][[eval(segcol)]] <- rev(lg[[i]][[eval(segcol)]])
      }
      # if user requested sections to be colored make list for this lg
      if (!is.null(sectcoldf))
      {
        sectcol[[i]] <- rev(subset(sectcoldf, sectcoldf$chr == lg[[i]][1, 1]))
      }
    }
    else {
      lg[[i]]$position <- round(lg[[i]]$position, roundpos)
      # if user requested sections to be colored make list for this lg
      if (!is.null(sectcoldf))
      {
        sectcol[[i]] <- subset(sectcoldf, sectcoldf$chr == lg[[i]][1, 1])
      }
    }
    
    # add to the list for this linkage group the qtls
    # and reverse postions if necessary
    
    if (!is.null(qtldf))
    {
      qtl[[i]] <- subset(qtldf, qtldf$chr == lg[[i]][1, 1])
      for (qtlix in 1:nrow(qtl[[i]])) {
        if (nrow(subset(qtldf, qtldf$chr == lg[[i]][1, 1])) != 0) {
          if (qtl[[i]]$chr[qtlix] %in% revthese) {
            # reverse positions
            so <-
              qtl[[i]]$so[qtlix]
            si <-
              qtl[[i]]$si[qtlix]
            ei <-
              qtl[[i]]$ei[qtlix]
            eo <-
              qtl[[i]]$eo[qtlix]
            qtl[[i]]$eo[qtlix] <-
              max(lg[[i]]$position) - so - min(lg[[i]]$position)
            qtl[[i]]$ei[qtlix] <-
              max(lg[[i]]$position) - si - min(lg[[i]]$position)
            qtl[[i]]$si[qtlix] <-
              max(lg[[i]]$position) - ei - min(lg[[i]]$position)
            qtl[[i]]$so[qtlix] <-
              max(lg[[i]]$position) - eo - min(lg[[i]]$position)
          }
        }
      }
    }
    
    # if user requested sections to be colored make list for this lg
    if (!is.null(sectcoldf))
    {
      sectcol[[i]] <- subset(sectcoldf, sectcoldf$chr == lg[[i]][1, 1])
    }
  }
  
  
  # -- Begin determine smallest and largest position for each row --------
  
  miny <- vector(length = nbrrows)
  maxy <- vector(length = nbrrows)
  
  
  for (nr in 1:nbrrows) {
    fromlg <- (nr - 1) * lgperrow + 1
    tolg <- min(c(length(mapthese),  nr * lgperrow))
    
    miny[nr] <- min(lg[[fromlg]]$position)
    maxy[nr] <- max(lg[[fromlg]]$position)
    
    for (i in fromlg:tolg)     {
      if (max(lg[[i]]$position) > maxy[nr]) {
        maxy[nr] <- max(lg[[i]]$position)
      }
      if (min(lg[[i]]$position) < miny[nr]) {
        miny[nr] <- min(lg[[i]]$position)
      }
    }
  }
  
  # -- End determine smallest and largest position for each row ----------
  
  # -- Begin loop for lg - apply formats & get required size -------------
  for (i in 1:length(mapthese)) {
    #if any of these options specified - apply first
    lrcol[[i]] <- rcol
    lrcol[[i]] = rep(lrcol[[i]], length.out = length(lg[[i]]$locus))
    
    lrfont[[i]] <- rfont
    lrfont[[i]] = rep(lrfont[[i]], length.out = length(lg[[i]]$locus))
    
    lrcex[[i]] <- rcex
    lrcex[[i]] <-
      rep(lrcex[[i]], length.out = length(lg[[i]]$locus))
    
    llcol[[i]] <- lcol
    llcol[[i]] = rep(llcol[[i]], length.out = length(lg[[i]]$position))
    
    llfont[[i]] <- lfont
    llfont[[i]] = rep(llfont[[i]], length.out = length(lg[[i]]$position))
    
    llcex[[i]] <- lcex
    llcex[[i]] <-
      rep(llcex[[i]], length.out = length(lg[[i]]$position))
    
    # next apply specific requst by locus
    if (!is.null(markerformatlist)) {
      for (ml in 1:length(markerformatlist)) {
        if (!is.null(markerformatlist[[ml]]$col) &&
            !is.na(markerformatlist[[ml]]$col)) {
          lrcol[[i]][which(lg[[i]]$locus %in% markerformatlist[[ml]]$locus)] <-
            markerformatlist[[ml]]$col
        }
        if (!is.null(markerformatlist[[ml]]$cex) &&
            !is.na(markerformatlist[[ml]]$cex)) {
          lrcex[[i]][which(lg[[i]]$locus %in% markerformatlist[[ml]]$locus)] <-
            markerformatlist[[ml]]$cex
        }
        if (!is.null(markerformatlist[[ml]]$font) &&
            !is.na(markerformatlist[[ml]]$font)) {
          lrfont[[i]][which(lg[[i]]$locus %in% markerformatlist[[ml]]$locus)] <-
            markerformatlist[[ml]]$font
        }
      }
    }
    if (!is.null(qtldf)) {
      qtldfone <- qtl[[i]]
      if (nrow(qtl[[i]]) == 0) {
        qtldfone <- NULL
      }
    }
    else {
      qtldfone <- NULL
    }
    
    
    dim[[i]] <- reqdim(
      lg[[i]],
      c(miny[nr], maxy[nr]),
      denmap = denmap,
      maxnbrcolsfordups = maxnbrcolsfordups,
      pdf.width = 30,
      labdist = labdist,
      lcol = llcol[[i]],
      lfont = llfont[[i]],
      lcex = llcex[[i]],
      rcol = lrcol[[i]],
      rfont = lrfont[[i]],
      rcex = lrcex[[i]],
      cex.lgtitle = cex.lgtitle,
      qtldf = qtldfone,
      ruler = ruler,
      prtlgtitles = prtlgtitles,
      showonly = showonly
    )
    
  }
  
  # -- End loop for lg - apply formats & get required size -------------
  
  pgxlg <- vector(length = length(mapthese))
  width <- vector(length = length(mapthese))
  
  totwidth <- rep(0, nbrrows)
  totheight <- rep(0, nbrrows)
  
  for (nr in 1:nbrrows) {
    fromlg <- (nr - 1) * lgperrow + 1
    tolg <- min(c(length(mapthese),  nr * lgperrow))
    
    #    calculate relative widths of plot areas
    #    and y range for all to be mapped
    
    # give a little for left so linkage group not next to ruler
    
    dim[[fromlg]]$reqwidth <-
      dim[[fromlg]]$reqwidth + 0.19
    #   dim[[tolg]]$reqwidth <-
    #    dim[[tolg]]$reqwidth + 0.3
    
    for (i in fromlg:tolg)    {
      totwidth[nr] <- dim[[i]]$reqwidth + totwidth[nr]
      if (dim[[i]]$reqheight > totheight[nr]) {
        totheight[nr] <- dim[[i]]$reqheight
      }
    }
  }
  
  allrowwidth <- max(totwidth)
  if (denmap & !is.null(sectcoldf$dens)) {
    # add a one 1/2 inch row for density map legend
    allrowheight <- sum(totheight) + 1.5
    relheight <- vector(length = nbrrows + 1)
    relheight[(nbrrows + 1)] <- 1.5/allrowheight
  }
  else {
    allrowheight <- sum(totheight)
    relheight <- vector(length = nbrrows)
  }
  
  # determine relative height of each row for layout
  
  
  for (nr in 1:nbrrows) {
    fromlg <- (nr - 1) * lgperrow + 1
    tolg <- min(c(length(mapthese),  nr * lgperrow))
    
    maxrowheight <- 0
    for (i in fromlg:tolg)    {
      if (dim[[i]]$reqheight > maxrowheight) {
        maxrowheight <- dim[[i]]$reqheight 
      }
    }
    relheight[nr] <- (maxrowheight / allrowheight) 
  }
  
  # add in margins
  allrowwidth <- allrowwidth + par("mai")[2] + par("mai")[4] + par("omi")[2] + par("omi")[4]
  allrowheight <- allrowheight + par("omi")[2] + par("omi")[4]
  
  # if user did not specify size, set to required size
  if (is.null(pdf.width)) {
    pdf.width <- ceiling(allrowwidth)
  }
  if (is.null(pdf.height)) {
    pdf.height <- ceiling(allrowheight)
  }
  
  message(c("Required pdf.width = ", allrowwidth))
  message(c("Required pdf.height = ", allrowheight))
  
  message(c("Using pdf.width = ", pdf.width))
  message(c("Using pdf.height = ", pdf.height))
  
  # turn off pdf used just for sizing and start the real one
  dev.off()
  pdf.options(
    bg = pdf.bg,
    title = pdf.title,
    family = pdf.family,
    pointsize = pdf.pointsize,
    fg = pdf.fg
  )
  
  pdf(outfile, width = pdf.width, height = pdf.height)
  
  par(mar = c(0, leftmar, topmar, 0),oma = c(1, 0.2, topoma, 0.2),new=FALSE)
  
  if (denmap & !is.null(sectcoldf$dens) ) {
    layout(c(seq(1, nbrrows + 1)), heights = relheight)
  }
  else {
    layout(c(seq(1, nbrrows)), heights = relheight)
  }
  
  # --- Begin loop for nbr rows to draw linkage groups -----------------
  for (nr in 1:nbrrows) {
    plot(
      .5,
      .5,
      xlim = c(0, 1),
      ylim = c(maxy[nr] + 7, miny[nr]),
      type = "n",
      cex = 1,
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      xaxs = "i",
      # don't pad x axis on each side
      bty = "n"
    )
    
    pin <- par("pin")[1]
    
    if (ruler) {
      axis(
        side = 2,
        col.axis = col.axis,
        cex.axis = cex.axis,
        font.axis = font.axis,
        at = at.axis,
        labels = labels.axis,
        lty = lty.axis,
        lwd = lwd.axis,
        lwd.ticks = lwd.ticks.axis
      )
      mtext(substitute(paste(bold(ylab))), side = 2, line=leftmar-0.8)
      # axis(
      #   side = 1,
      #   col.axis = col.axis,
      #   cex.axis = cex.axis,
      #   font.axis = font.axis,
      #   at = seq(0.05, 1.05, 0.066667),
      #   labels = c("1", "2", "3", "4", "5", 
      #              "6", "7", "8", "9", "10", 
      #              "11", "12", "13", "14", "15"),
      #   lty = lty.axis,
      #   lwd = lwd.axis,
      #   outer = T,
      #   pos = 120,
      #   lwd.ticks = lwd.ticks.axis
      # )
      # mtext("Linkage group", side = 1, line=1)
      
    }
    
    # if last row put footnote
    # if (nr == nbrrows) {
    #   mtext(
    #     "Rendered by LinkageMapView",
    #     side = 1,
    #     cex = 0.5,
    #     outer=TRUE,
    #     col = pdf.fg,
    #     adj = 0
    #   )
    # }
    # set up list for returned values from drawone
    # needed for placement of segments connecting markers
    # between linkage groups
    
    yrlabwidth <- list()
    adjyr <- list()
    adjyl <- list()
    dups <- list()
    
    
    widthused <- 0
    
    fromlg <- (nr - 1) * lgperrow + 1
    tolg <- min(c(length(mapthese),  nr * lgperrow))
    
    # ----------- Begin loop for each lg in this row -----------------------
    for (i in fromlg:tolg) {
      pctwidth <- dim[[i]]$reqwidth / totwidth[nr]
      width <- pctwidth * pin
      if (posonleft[i]) {
        if (!ruler) {
          llablen <- dim[[i]]$maxlenllab
          lspace <-
            strwidth("M", units = "inches") * max(lcex)
          llabdist <- labdist
        }
        else {
          if (i > 1 && !posonleft[i - 1]) {
            llablen <- 0
            # make sure titles do not overlap when mirror no labels come together
            lspace <-
              max(strwidth(lg[[i]][1, 1]) * cex.lgtitle ,
                  strwidth(lg[[i - 1]][1, 1]) * cex.lgtitle) / 2 + strwidth("M", units = "inches") * cex.lgtitle
            llabdist <- 0
          }
          else {
            llablen <- 0
            lspace <- 0
            llabdist <- 0
          }
        }
      }
      else {
        llablen <- dim[[i]]$maxlenrlab
        lspace <-
          strwidth("M", units = "inches") * max(rcex)
        llabdist <- labdist
      }
      pgxlg[i] <-
        (sum(llablen,
             lgw / 2,
             llabdist,
             lspace)) / pin  + widthused / pin
      
      # give a little for left margin on first one on row
      if (i == (nr-1)*lgperrow +1) {
        pgxlg[i] <- pgxlg[i] + 0.3 / pdf.width
      }
      widthused <- widthused + width
      
      if (!is.null(qtldf)) {
        qtldfone <- qtl[[i]]
        if (nrow(qtl[[i]]) == 0) {
          qtldfone <- NULL
        }
      }
      else {
        qtldfone <- NULL
        
      }
      
      if (!is.null(sectcoldf)) {
        sectcoldfone <- sectcol[[i]]
        if (nrow(sectcol[[i]]) == 0) {
          sectcoldfone <- NULL
        }
      }
      else {
        sectcoldfone <- NULL
      }
      
      
      if (!is.null(lgtitles)) {
        lgtitleone <- lgtitles[i]
      }
      else {
        lgtitleone <- NULL
      }
      
      if (i == 1) {
        #pass title to drawone
        main <- main
      }
      else {
        main = NULL
      }
      
      dolist <-
        drawone(
          lg[[i]],
          dim[[i]],
          totwidth[nr],
          c(miny[nr], maxy[nr]),
          denmap = denmap,
          maxnbrcolsfordups = maxnbrcolsfordups,
          pdf.width = pin,
          pdf.fg = pdf.fg,
          lgw = lgw,
          lg.col = lg.col,
          lg.lwd = lg.lwd,
          pgx = pgxlg[i] ,
          labdist = labdist ,
          lcol = llcol[[i]],
          lfont = llfont[[i]],
          lcex = llcex[[i]],
          rcol = lrcol[[i]],
          rfont = lrfont[[i]],
          rcex = lrcex[[i]],
          rsegcol = rsegcol,
          main = main,
          cex.main = cex.main,
          font.main = font.main,
          col.main = col.main,
          cex.lgtitle = cex.lgtitle,
          font.lgtitle = font.lgtitle,
          col.lgtitle = col.lgtitle,
          qtldf = qtldfone,
          posonleft = posonleft[i],
          ruler = ruler,
          prtlgtitles = prtlgtitles,
          lgtitles = lgtitleone,
          segcol = segcol,
          showonly = showonly,
          sectcoldf = sectcoldfone
        )
      yrlabwidth[[i]] <- dolist$yrlabwidth
      adjyr[[i]] <- dolist$adjyr
      adjyl[[i]] <- dolist$adjyl
      dups[[i]] <- dolist$dups
      solist[[i]] <- dolist$solist
      
    }
    
    # ----------- End loop for each lg in this row -----------------------
    
    # connect markers across linkage groups if requested
    
    # make sure connect marker data frame passed has no factors
    
    if (autoconnadj == TRUE)
    {
      #only pass the markers on this row
      
      autoconndf <- autoconn(lg[fromlg:tolg], pdf.fg, lgperrow)
      if (!is.null(conndf)) {
        fas <- sapply(conndf, is.factor)
        conndf[fas] <- lapply(conndf[fas], as.character)
        allconndf <- rbind(conndf, autoconndf)
        # get rid of duplicates if user specified and automatically
        # since auto is at the end, the user specified will be kept
        allconndf <- allconndf[!duplicated(allconndf[, 1:4]), ]
      }
      else {
        allconndf <- autoconndf
      }
    }
    else {
      if (!is.null(conndf)) {
        fas <- sapply(conndf, is.factor)
        conndf[fas] <- lapply(conndf[fas], as.character)
        allconndf <- conndf
      }
      else{
        allconndf <- data.frame() # empty
      }
    }
    
    if (nrow(allconndf) > 0) {
      for (i in 1:nrow(allconndf)) {
        # determine index for from and to chr
        
        fromi <- match(allconndf$fromchr[i], mapthese)
        if (is.na(fromi)) {
          stop(
            c(
              "Connect marker from position not found for chr = ",
              allconndf$fromchr[i],
              " and locus = ",
              allconndf$fromlocus[i]
            )
          )
        }
        toi <- match(allconndf$tochr[i], mapthese)
        if (is.na(toi)) {
          stop(
            c(
              "Connect marker to position not found for chr = ",
              allconndf$tochr[i],
              " and locus = ",
              allconndf$tolocus[i]
            )
          )
        }
        if (toi < fromi) {
          temp <- toi
          toi <- fromi
          fromi <- temp
        }
        
        
        # user requested connections must be in same row
        if (ceiling(fromi / lgperrow) != ceiling(toi / lgperrow)) {
          stop(
            c(
              "Connect marker from chr ",
              allconndf$fromchr[i],
              " must be on the same row as to chr ",
              allconndf$tochr[i],
              " and you have lgperrow = ",
              lgperrow
            )
          )
        }
        
        # determine x and y for from and to marker
        # look up marker in lgin to get position
        fpos <-
          lg[[fromi]]$position[lg[[fromi]]$locus == allconndf$fromlocus[i]]
        if (identical(fpos, numeric(0))) {
          stop(
            c(
              "Connect marker from position not found for chr = ",
              allconndf$fromchr[i],
              " and locus = ",
              allconndf$fromlocus[i]
            )
          )
        }
        tpos <-
          lg[[toi]]$position[lg[[toi]]$locus == allconndf$tolocus[i]]
        if (identical(tpos, numeric(0))) {
          stop(
            c(
              "Connect marker to position not found for chr = ",
              allconndf$tochr[i],
              " and locus = ",
              allconndf$tolocus[i]
            )
          )
        }
        
        if (!is.null(showonly)) {
          fy <- match(fpos, solist[[fromi]]$newllab)
          ty <- match(tpos, solist[[toi]]$newllab)
        }
        else
        {
          fy <- match(fpos, lg[[fromi]]$position)
          ty <- match(tpos, lg[[toi]]$position)
        }
        # determine from
        if (posonleft[fromi]) {
          connfxpos <-
            ((yrlabwidth[[fromi]][fy] + lgw / 2 + labdist / 2) / pin) + pgxlg[fromi]
          if (!is.null(showonly)) {
            connfypos <-
              adjyr[[fromi]][match(fpos, solist[[fromi]]$newllab[dups[[fromi]]$rkeep])]
          }
          else {
            connfypos <-
              adjyr[[fromi]][match(fpos, lg[[fromi]]$position[dups[[fromi]]$rkeep])]
          }
        }
        
        else {
          if (!ruler) {
            connfxpos <-
              ((
                dim[[fromi]]$maxlenllab + lgw / 2 + labdist + strwidth("M", units = "inches")
              ) / pin) + pgxlg[fromi]
            if (!is.null(showonly)) {
              connfypos <-
                adjyl[[fromi]][match(fpos, solist[[fromi]]$newllab[dups[[fromi]]$lkeep])]
            }
            else {
              connfypos <-
                adjyl[[fromi]][match(fpos, lg[[fromi]]$position[dups[[fromi]]$lkeep])]
            }
            
          } else
          {
            connfxpos <- ((lgw / 2) / pin) + pgxlg[fromi]
            connfypos <-
              fpos
            
          }
        }
        
        if (!posonleft[toi]) {
          conntxpos <-
            -((yrlabwidth[[toi]][ty]  + lgw / 2 + labdist / 2) / pin) + pgxlg[toi]
          if (!is.null(showonly)) {
            conntypos <-
              adjyr[[toi]][match(tpos, solist[[toi]]$newllab[dups[[toi]]$rkeep])]
          }
          else {
            conntypos <-
              adjyr[[toi]][match(tpos, lg[[toi]]$position[dups[[toi]]$rkeep])]
          }
        }
        
        else {
          if (!ruler) {
            conntxpos <-
              -((
                dim[[toi]]$maxlenllab  +  lgw / 2 + labdist + strwidth("M", units = "inches")
              ) / pin) + pgxlg[toi]
            if (!is.null(showonly)) {
              conntypos <-
                adjyl[[toi]][match(tpos, solist[[toi]]$newllab[dups[[toi]]$lkeep])]
            }
            else {
              conntypos <-
                adjyl[[toi]][match(tpos, lg[[toi]]$position[dups[[toi]]$lkeep])]
            }
          }
          else {
            conntxpos <- -((lgw / 2) / pin) + pgxlg[toi]
            conntypos <-
              tpos
            
          }
        }
        if (is.null(allconndf$col[i])) {
          allconndf$col[i] = pdf.fg
        }
        segments(connfxpos,
                 connfypos,
                 conntxpos,
                 conntypos,
                 col = allconndf$col[i])
      }
    }
  }
  
  
  # --- End loop for nbr rows to draw linkage groups -----------------
  
  # if density map legend to be displayed
  if (denmap &
      !is.null(sectcoldf$dens)) {
    # plot legend is last row
    # na.last removes NAs
    leg <- sectcoldf[order(sectcoldf$dens,na.last=NA), ]
    if (max(leg$dens) > 100) {
      leg$dens <- round(leg$dens, digits = 0)
    }
    else {
      leg$dens <- round(leg$dens, digits = 1)
    }
    uleg <- leg[!duplicated(leg$dens), ]
    # reduce to reasonable number of buckets 1/4 inch wide bucket minimum
    if ((pdf.width / .25) < nrow(uleg)) {
      nbrbuckets <- round(pdf.width / .25, digits = 0)
      bplotdens <-
        uleg$dens[seq(1, length(uleg$dens), length(uleg$dens) / nbrbuckets)]
      bplotcol <-
        uleg$col[seq(1, length(uleg$col), length(uleg$col) / nbrbuckets)]
    } else {
      bplotdens <- uleg$dens
      bplotcol <- uleg$col
    }
    if (!bplotdens[length(bplotdens)] == uleg$dens[length(uleg$dens)]) {
      # include largest density
      bplotdens <- append(bplotdens, uleg$dens[length(uleg$dens)])
      bplotcol <- append(bplotcol, uleg$col[length(uleg$col)])
    }
    par(mar = c(5, leftmar, 1, 1))
    barplot(
      rep(0.1, length(bplotcol)),
      col = bplotcol,
      space = 0,
      axes = F,
      xlab = paste("Marker interval (", units, "/Locus)", sep = ""),
      names = bplotdens,
      cex.names = .75,
      las = 2,
      cex.lab = .75
    )
    
  }
}

#############################################################################################
#############################################################################################
convertfont  <- function (parm, fonttext) {
  # make sure valid
  validfont <-
    c("plain text", "bold", "italic", "bold italic", "1", "2", "3", "4")
  retfont <- vector(length = length(fonttext))
  
  for (i in 1:length(fonttext)) {
    if (!(fonttext[i] %in% validfont)) {
      message("Invalid font entered for paramenter ",
              parm, ":", fonttext[i])
      if (length(fonttext) > 1) {
        message ("at row number ", i)
      }
      message("Please use one of the following standard R fonts or see pdf.family parameter: ")
      message("1 or 'plain text'")
      message("2 or 'bold'")
      message("3 or 'italic'")
      message("4 or 'bold italic'")
      stop("")
    }
    else {
      if (fonttext[i] == "1" | fonttext[i] == "plain text") {
        retfont[i] <- 1
      }
      
      
      else {
        if (fonttext[i] == "2" | fonttext[i] == "bold") {
          retfont[i] <-  2
        }
        
        
        else {
          if (fonttext[i] == "3" | fonttext[i] == "italic") {
            retfont[i] <- 3
          }
          
          
          else {
            if (fonttext[i] == "4" | fonttext[i] == "bold italic") {
              retfont[i] <- 4
            }
          }
        }
      }
    }
  }
  
  retfont
}

readlgdf <-
  function(df,
           mapthese)
  {
    if (ncol(df) < 3) {
      for (i in 1:ncol(df)) {
        stop("Less than 3 columns found on input data frame - see help for mapthis parameter")
      }
      
    }
    
    # assign my own column names so I can reference by name
    colnames(df)[1:3] <- c("group", "position", "locus")
    
    # make group a character field
    
    df$group <- as.character(df$group)
    
    if (!missing(mapthese)) {
      if (!all(mapthese %in% df$group)) {
        stop ("chrnames to map not found in input data frame")
      }
    }
    
    # make sure data frame passed has no factors
    if (!is.null(df)) {
      fas <- sapply(df, is.factor)
      df[fas] <- lapply(df[fas], as.character)
    }
    
    return (df)
  }

getlg <- function(df, lg, dupnbr,roundpos) {
  # if dupnbr is true user wants only first marker name
  # at a position to show with (### more) after the marker name
  # to indicate how many more markers at that position
  # We'll use (### more) as the second label so that the first
  # label isn't corrupted for connecting homologous markers
  # based on matching marker names.  That means we need to
  # set maxnbrcolsfordups <- 2 in lmv.linkage.plot.
  
  templg <- subset(df, df[, 1] == lg)
  # make sure linkage group is in order by position
  thislg <- templg[order(templg$position),]
  thislg$position <- round(thislg$position, roundpos)
  if (dupnbr) {
    keeprows <- vector()
    dupcount <- vector(length = nrow(thislg))
    keeprows <- append(keeprows, 1)  # always keep first row
    dupcount[1] <- 0
    i <- 0
    for (i in 2:nrow(thislg)) {
      if (thislg$position[i] == thislg$position[i - 1]) {
        dupcount[i] <- dupcount[i - 1] + 1
        if (dupcount[i] == 1) {
          keeprows <- append(keeprows, i) # for (## more label)
        }
      }
      else {
        if (dupcount[i - 1] > 0) {
          thislg$locus[i - dupcount[i - 1]] <-
            paste(" (",
                  dupcount[i - 1],
                  " more)",
                  sep = "")
        }
        keeprows <- append(keeprows, i)
        dupcount[i] <- 0
      }
    }
    # finish last row
    if (dupcount[i] > 0) {
      thislg$locus[i - dupcount[i - 1]] <-
        paste(" (",
              dupcount[i],
              " more)",
              sep = "")
    }
    retdf <- thislg[unique(keeprows), ]
  }
  else {
    retdf <- thislg
  }
  return(retdf)
}

editlgdf  <- function (df) {
  # make sure all positions are numeric
  if (!is.numeric(df[,2]))
  { for (i in 1:nrow(df)) {
    message(c("Some positions are not numeric - in linkage group: ",df[i,1]))
    if (!is.numeric(df[i,2])) {
      message(c("I see ",  df[i,2] , " in row " , i))
      
    }
  }
    
    stop("Some positions are not numeric") }
}

reqdim <- function(df,
                   yrange,
                   pdf.width = 12,
                   denmap = FALSE,
                   maxnbrcolsfordups = 3,
                   lgw = 0.25,
                   pgx = 0.5 ,
                   labdist = 0.1 ,
                   rcex = par("cex"),
                   lcex = par("cex"),
                   rfont = par("font"),
                   lfont = par("font"),
                   rcol = par("col"),
                   lcol = par("col"),
                   cex.lgtitle = cex.lgtitle,
                   qtldf = NULL,
                   ruler = FALSE,
                   prtlgtitles = TRUE,
                   showonly = showonly)
{
  
  # set up points to label - the points will be plotted invisibly
  # x1 is on the left side of the chromosome, x2 on the right
  # these will be recalculated after required dimensions are determined
  
  lgwpct <- lgw / pdf.width
  x1 <-
    rep(pgx - lgwpct / 2, length(df$position))
  x2 <-
    rep(pgx + lgwpct / 2, length(df$position))
  
  y <- df$position
  rlab <- df$locus
  llab <- df$position
  
  # ylim reversed so 0 is at top of y-axis
  # don't print the points: type="n"
  # and don't print axis: xaxt and yaxt = "n"
  # and don't print axis labels: yaxt="n",xlab=""
  # and don't print the box around the plot: bty="n"
  # in other words don't print anything but establish the points
  # for the markers
  
  
  plot(
    x2,
    y,
    xlim = c(0, 1),
    ylim = rev(range(y)),
    type = "n",
    cex = 1,
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    xaxs = "i",
    bty = "n"
    
  )
  
  if (!is.null(showonly)) {
    solist <- show(
      showonly,
      llab,
      rlab,
      rcex = rcex,
      lcex = lcex,
      rfont = rfont,
      lfont = lfont,
      rcol = rcol,
      lcol = lcol
    )
    llab <- solist$newllab
    rlab <- solist$newrlab
    rcex <- solist$newrcex
    lcex <- solist$newlcex
    rfont <- solist$newrfont
    lfont <- solist$newlfont
    rcol <- solist$newrcol
    lcol <- solist$newlcol
  }
  
  yrlabwidth <- vector(length = length(llab))
  if (!denmap) {
    # find dups to figure out how many columns to reserve space for
    dups <- fsdups(llab, maxnbrcolsfordups)
    
    # Determine width of each right label including dups
    
    yrlabwidth[setdiff(dups$rkeep, dups$frkeep)] <-
      strwidth(rlab[setdiff(dups$rkeep, dups$frkeep)], units = "inches") *
      rcex[setdiff(dups$rkeep, dups$frkeep)] + labdist
    yrlabwidth[dups$frkeep] <-
      strwidth(rlab[dups$frkeep], units = "inches") *
      rcex[dups$frkeep] + labdist * 1.2
    # now add space for duplicates in cols
    if (maxnbrcolsfordups > 1) {
      for (i in 1:length(llab)) {
        for (m in 1:(maxnbrcolsfordups - 1)) {
          if (!is.na(dups$yd[m, i])) {
            yrlabwidth[i] <-
              yrlabwidth[i] + strwidth(rlab[(i + m)], units = "inches") * rcex[(i + m)] + (strwidth(" ", units =
                                                                                                      "inches")  * rcex[(i + m)]) / 2
          }
        }
      }
    }
    
    # If qtl provided, calculate width required
    # Since we don't know pdf dimensions yet, we
    # can't determine if qtl will overlap each other.
    # so we assume they will and add enough space to
    # the required width for each qtl
    
    # adding in strwidth("M",units=inches) takes care of
    # 0.5 of a character width default between points and labels
    # and leaves a little room between linkage groups and at edge of pdf
    if (!ruler) {
      reqwidth <- sum(
        max(strwidth(llab, units = "inches") * lcex),
        max(yrlabwidth),
        lgw,
        labdist * 2,
        strwidth("M", units = "inches")  * max(lcex),
        strwidth("M", units  = "inches")  * max(rcex),
        nrow(qtldf) * (lgw / 3 + strheight("M", units = "inches") *
                         3)
      )
      
      reqheight <-
        max(
          sum(strheight(llab[dups$rkeep], units = "inches") * lcex[dups$rkeep] * 1.4),
          sum(strheight(rlab[dups$rkeep], units = "inches") * rcex[dups$rkeep] * 1.4),
          sum(strheight(llab[dups$rkeep], units = "inches") * rcex[dups$rkeep] * 1.4)
        ) # because positions spread like markers
    }
    else {
      reqwidth <- sum(
        max(yrlabwidth),
        lgw ,
        labdist,
        strwidth("M", units  = "inches")  * max(rcex),
        nrow(qtldf) * (lgw / 3 + strheight("M", units = "inches") *
                         3)
      )
      
      reqheight = sum(strheight(rlab, units = "inches") * rcex * 1.4)
      
    }
  }
  # else we are doing a density map
  
  else {
    yrlabwidth <- strwidth("M", units = "inches")
    # only doing density map
    # line segments are 1/96"
    # so make sure closest two lines only half overlap
    reqheight <- (diff(range(y) / (2 * 96)) / min(diff(y)[diff(y) > 0]))
    # width is max of title or chromosome plus any qtls
    reqwidth <-
      sum(
        max(
          lgw ,
          strwidth(df$group, units = "inches") * cex.lgtitle + strwidth("M", units = "inches") *
            cex.lgtitle
        ),
        nrow(qtldf) * (lgw / 3 + strheight("M", units = "inches") *
                         3)
      )
  }
  
  # give a margin at top and bottom for chromosome ends and margins
  reqheight = reqheight + lgw + par("mai")[1] + par("mai")[3]
  
  list(
    reqwidth = reqwidth,
    reqheight = reqheight,
    maxlenrlab = max(yrlabwidth),
    maxlenllab = max(strwidth(llab, units = "inches") * lcex)
  )
  
}

autoconn <-
  function(lg, fg, lgperrow) {
    conndf <- data.frame(
      fromchr = character(),
      fromlocus = character(),
      tochr = character(),
      tolocus = character(),
      col = character(),
      stringsAsFactors = FALSE
    )
    
    # first remove any "(## more)" labels generated by dupnbr=TRUE
    newlg <- list()
    for (i in 1:length(lg)) {
      morelab <- lg[[i]]$locus[grep("more)", lg[[i]]$locus)]
      newlg[[i]] <-
        subset(lg[[i]],!(lg[[i]]$locus %in% morelab) &
                 !(lg[[i]]$locus == ""))
    }
    
    # next pick out the from
    if (length(newlg) > 1) {
      for (i in 2:length(newlg)) {
        # if no locus labels at all skip it (density map)
        if (nrow(newlg[[i - 1]]) > 0 &
            nrow(newlg[[i]]) > 0) {
          # don't connect if for some reason they requested
          # the same linkage group side by side
          # or if they are going to be in different rows
          if (!(newlg[[i - 1]][1, 1] == newlg[[i]][1, 1])
              && (i - 1) %% lgperrow != 0) {
            connlocus <-
              newlg[[i - 1]]$locus[newlg[[i - 1]]$locus %in% newlg[[i]]$locus]
            if (length(connlocus) > 0) {
              conndf <- rbind(
                conndf,
                data.frame(
                  fromchr = newlg[[i - 1]][1, 1],
                  fromlocus = connlocus,
                  tochr = newlg[[i]][1, 1],
                  tolocus = connlocus,
                  col = fg
                )
              )
            }
          }
        }
      }
    }
    return (conndf)
  }

drawone <-
  function(df,
           dim,
           totwidth,
           yrange,
           denmap = FALSE,
           maxnbrcolsfordups = 3,
           pdf.width = 12,
           pdf.fg = "black",
           lgw = 0.25,
           lg.col = NULL,
           lg.lwd = par("lwd"),
           pgx = 0.5 ,
           labdist = 0.1 ,
           rcex = par("cex"),
           lcex = par("cex"),
           rfont = par("font"),
           lfont = par("font"),
           rcol = par("col"),
           lcol = par("col"),
           rsegcol = TRUE,
           main = main,
           cex.main = cex.main,
           font.main = font.main,
           col.main = col.main,
           cex.lgtitle = par("cex.main"),
           font.lgtitle = par("font.main"),
           col.lgtitle = par("col.main"),
           qtldf = NULL,
           posonleft = TRUE,
           ruler = FALSE,
           prtlgtitles = TRUE,
           lgtitles = NULL,
           segcol = NULL,
           showonly = NULL,
           sectcoldf = NULL) {
    y <- df$position
    rlab <- df$locus
    llab <- df$position
    
    pctwidth <- dim$reqwidth / totwidth
    width <- pctwidth * pdf.width
    
    lgwpct <- lgw / pdf.width
    
    # if user requested to have positions show on right
    # set up adjustments
    
    
    if (posonleft) {
      labdistpct <- labdist / pdf.width
      x1 <- rep(pgx - lgwpct / 2, length(df$position))
      x2 <- rep(pgx + lgwpct / 2, length(df$position))
      rpos <- 4
      lpos <- 2
      posmult <- 1    # adjuster for left or right positioning
    }
    else {
      labdistpct <- -labdist / pdf.width
      x1 <- rep(pgx + lgwpct / 2, length(df$position))
      x2 <- rep(pgx - lgwpct / 2, length(df$position))
      rpos <- 2
      lpos <- 4
      posmult <- -1
    }
    
    if (!is.null(segcol)) {
      linesegcolor <- df[[eval(segcol)]]
    }
    else {
      linesegcolor <- rep(pdf.fg,nrow(df))
    }
    
    points(
      x2,
      y,
      type = "n",
      cex = 1,
      xlab = "",
      ylab = ""
    )
    
    # put title on top of chromosome
    if (prtlgtitles) {
      if (!is.null(lgtitles)) {
        lgtext = lgtitles
      }
      else {
        lgtext = df[1, 1]
      }
      #mtext writes text in the margins
      
      mtext(
        lgtext,
        at = pgx,
        line = 0.05,
        cex = cex.lgtitle,
        col = col.lgtitle,
        font = font.lgtitle
      )
    }
    if (!is.null(main)) {
      mtext(
        main,
        at = .5,
        line = 1,
        outer = TRUE,
        cex = cex.main,
        col = col.main,
        font = font.main
      )
    }
    
    
    # eliminate all but showonly labels if requested
    
    if (!is.null(showonly)) {
      solist <- show(
        showonly,
        llab,
        rlab,
        rcex = rcex,
        lcex = lcex,
        rfont = rfont,
        lfont = lfont,
        rcol = rcol,
        lcol = lcol
      )
      llab <- solist$newllab
      rlab <- solist$newrlab
      rcex <- solist$newrcex
      lcex <- solist$newlcex
      rfont <- solist$newrfont
      lfont <- solist$newlfont
      rcol <- solist$newrcol
      lcol <- solist$newlcol
    }
    else {
      solist <- NULL
    }
    # if density map skip all of this
    if (!denmap) {
      # find and save dup locations before calling spreadcexlabs
      dups <- fsdups(llab, maxnbrcolsfordups)
      
      # spread the labels as necessary
      # heights are negative because y axis is reversed
      
      if (length(dups$rkeep) > 1) {
        adjyr <- spreadcexlabs(
          llab[dups$rkeep],
          max(strheight(rlab)) * -1.4,
          strh = -max(strheight(rlab)),
          min = min(yrange),
          max = max(yrange),
          cex = rcex[dups$rkeep],
          maxiter = 999999
        )
        
      }
      else {
        adjyr <- llab[dups$rkeep]
      }
      
      pos = rpos
      
      # label the first columns except those that are dups(
      if (length(setdiff(dups$rkeep, dups$frkeep)) > 0) {
        text(
          x2[setdiff(dups$rkeep, dups$frkeep)] + labdistpct,
          adjyr[setdiff(1:length(adjyr), dups$fykeep)],
          labels = rlab[setdiff(dups$rkeep, dups$frkeep)],
          pos = pos,
          col = rcol[setdiff(dups$rkeep, dups$frkeep)],
          cex = rcex[setdiff(dups$rkeep, dups$frkeep)],
          font = rfont[setdiff(dups$rkeep, dups$frkeep)]
        )
      }
      
      # label the first columns of dups inset for visibility
      if (length(dups$frkeep > 0)) {
        text(
          x2[dups$frkeep] + labdistpct * 1.2,
          adjyr[dups$fykeep],
          labels = rlab[dups$frkeep],
          pos = pos,
          col = rcol[dups$frkeep],
          cex = rcex[dups$frkeep],
          font = rfont[dups$frkeep]
        )
      }
      
      # now fill out duplicates in columns
      if (maxnbrcolsfordups > 1) {
        for (i in 1:length(llab)) {
          for (m in 1:(maxnbrcolsfordups - 1)) {
            if (!is.na(dups$yd[m, i])) {
              if (m == 1) {
                rx <- x2[i] + labdistpct * 1.2
              }
              ry <- adjyr[dups$yd[m, i]]
              rx <-
                rx + posmult * (strwidth(rlab[(i + m - 1)]) * rcex[(i + m - 1)] + (strwidth(" ")  * rcex[(i + m -
                                                                                                            1)]) / 2)
              
              text(
                rx,
                ry,
                labels = rlab[(i + m)],
                pos = pos,
                col = rcol[(i + m)],
                cex = rcex[(i + m)],
                font =
                  rfont[(i + m)]
              )
            }
          }
        }
      }
    }  # end skip all of this if denmap
    
    # xpd = NA to turn off clipping and arcs can go into margins
    par(xpd = NA)
    if (!is.null(lg.col)) {
      rect(pgx - lgwpct / 2, min(y), pgx + lgwpct / 2, max(y), col = lg.col)
      symbols(
        x = pgx,
        y = min(y),
        circles = lgwpct / 2,
        bg = lg.col,
        add = TRUE,
        fg = lg.col,
        inches = FALSE
      )
      symbols(
        x = pgx,
        y = max(y),
        circles = lgwpct / 2,
        bg = lg.col,
        add = TRUE,
        fg = lg.col,
        inches = FALSE
      )
    }
    
    
    # color sections and color arcs at end same as first and last section
    if (!is.null(sectcoldf)) {
      if (nrow(sectcoldf) > 0) {
        if (is.null(lg.col)) {
          if (denmap)
            #lg.col overrides coloring same as adjcent color
          {
            symbols(
              x = pgx,
              y = min(y),
              circles = lgwpct / 2,
              bg = sectcoldf$col[1],
              add = TRUE,
              fg = sectcoldf$col[1],
              inches = FALSE
            )
            symbols(
              x = pgx,
              y = max(y),
              circles = lgwpct / 2,
              bg = sectcoldf$col[nrow(sectcoldf)],
              add = TRUE,
              fg = sectcoldf$col[nrow(sectcoldf)],
              inches = FALSE
            )
          }
          for (sc in 1:nrow(sectcoldf)) {
            rect(
              pgx - lgwpct / 2,
              sectcoldf$s,
              pgx + lgwpct / 2,
              sectcoldf$e,
              col = sectcoldf$col,
              border = NA
            )
          }
        }
      }
    }
    
    segments(pgx - lgwpct / 2, min(y), pgx - lgwpct / 2, max(y), lwd = lg.lwd)
    segments(pgx + lgwpct / 2, min(y), pgx + lgwpct / 2, max(y), lwd = lg.lwd)
    
    plotrix::draw.arc(
      x = pgx,
      y = min(y),
      radius = lgwpct / 2,
      deg1 = 0,
      deg2 = 180,
      lwd = lg.lwd
    )
    plotrix::draw.arc(
      x = pgx,
      y = max(y),
      radius = lgwpct / 2,
      deg1 = -180,
      deg2 = 0,
      lwd = lg.lwd
    )
    
    if (!denmap) {
      if (rsegcol) {
        segcolprt <- rcol[setdiff(dups$rkeep, dups$frkeep)]
      }
      else {
        segcolprt <- pdf.fg
      }
      # segments for nondups from chr to marker
      
      segments(x2[setdiff(dups$rkeep, dups$frkeep)] + labdistpct,
               adjyr[setdiff(1:length(adjyr), dups$fykeep)],
               x2[setdiff(dups$rkeep, dups$frkeep)],
               llab[setdiff(dups$rkeep, dups$frkeep)],
               col = segcolprt)
      
      if (rsegcol) {
        segcolprt <- rcol[setdiff(dups$rkeep, dups$frkeep)]
      }
      else {
        segcolprt <- linesegcolor[setdiff(dups$rkeep, dups$frkeep)]
      }
      #connect across chromosome
      segments(x1[setdiff(dups$rkeep, dups$frkeep)],
               llab[setdiff(dups$rkeep, dups$frkeep)],
               x2[setdiff(dups$rkeep, dups$frkeep)],
               llab[setdiff(dups$rkeep, dups$frkeep)], col = segcolprt)
      
      if (rsegcol) {
        segcolprt <- rcol[dups$frkeep]
      }
      else {
        segcolprt <- pdf.fg
      }
      #segments for dups
      if (length(dups$frkeep) > 0) {
        segments(x2[dups$fykeep] + labdistpct * 1.2, adjyr[dups$fykeep],
                 x2[dups$fykeep],
                 llab[dups$frkeep],
                 col = segcolprt)
        
        if (rsegcol) {
          segcolprt <- rcol[dups$frkeep]
        }
        else {
          segcolprt <- linesegcolor[dups$frkeep]
        }
        
        #connect across chromosome
        segments(x1[dups$frkeep],
                 llab[dups$frkeep],
                 x2[dups$frkeep],
                 llab[dups$frkeep], col = segcolprt)
        
        #segments connecting dups vertically
        if (length(y[dups$rkeep]) > 1) {
          for (i in 2:length(y[dups$rkeep])) {
            if (llab[dups$rkeep][i] == llab[dups$rkeep][i - 1]) {
              segments(x2[dups$rkeep][i] + labdistpct * 1.2,
                       adjyr[i],
                       x2[dups$rkeep][i] + labdistpct * 1.2,
                       adjyr[(i - 1)])
            }
          }
        }
      }
      
      # now draw lines across chromosome for non-shown markers if any
      
      if (!is.null(showonly)) {
        segments(
          rep(x1[1], length.out = length(setdiff(y, (
            llab[dups$rkeep]
          )))),
          setdiff(y, (llab[dups$rkeep])),
          rep(x2[1], length.out = length(setdiff(y, (
            llab[dups$rkeep]
          )))),
          setdiff(y, (llab[dups$rkeep])),
          col = linesegcolor[match(setdiff(y, (llab[dups$rkeep])), df$position)]
        )
      }
      
      pos = lpos
      if (!ruler) {
        text(
          x1[dups$lkeep] - labdistpct,
          adjyr[dups$flkeep],
          labels = llab[dups$lkeep],
          pos = pos,
          col = lcol[dups$lkeep],
          cex = lcex[dups$lkeep],
          font =
            lfont[dups$lkeep]
        )
        
        segments(x1[dups$lkeep] - labdistpct, adjyr[dups$flkeep], x1[dups$lkeep], llab[dups$lkeep])
      }
      
    } # end don't do this for density map
    else {
      segments(rep(x1[1], length.out = length(y)),
               y,
               rep(x2[1], length.out = length(y)),
               y,
               col = linesegcolor)
    }
    
    # figure locus labwidth to pass back for connecting markers
    # and for drawing qtls
    
    
    if (denmap) {
      yrlabwidth <- rep(strwidth("M", units = "inches"),length(y))
      adjyr <- y
      adjyl <- y
      dups <- NULL
    }
    else {
      yrlabwidth <- vector(length = length(llab))
      yrlabwidth[setdiff(dups$rkeep, dups$frkeep)] <-
        strwidth(rlab[setdiff(dups$rkeep, dups$frkeep)], units = "inches") *
        rcex[setdiff(dups$rkeep, dups$frkeep)] + labdist
      yrlabwidth[dups$frkeep] <-
        strwidth(rlab[dups$frkeep], units = "inches") *
        rcex[dups$frkeep] + labdist * 1.2
      # now add space for duplicates in cols
      if (maxnbrcolsfordups > 1) {
        for (i in 1:length(llab)) {
          for (m in 1:(maxnbrcolsfordups - 1)) {
            if (!is.na(dups$yd[m, i])) {
              yrlabwidth[i] <-
                yrlabwidth[i] + strwidth(rlab[(i + m)], units = "inches") * rcex[(i + m)] + (strwidth(" ", units =
                                                                                                        "inches")  * rcex[(i + m)]) / 2
            }
          }
        }
      }
    }
    # save yrlabwidth before adding in any qtls to pass back for connecting markers
    returnlabwidth <- yrlabwidth
    
    # draw qtls
    if (!is.null(qtldf)) {
      # save calculated start and end (biggest of actual positions or label length)
      # for determining if labels will overlap
      sdf <- vector(length = nrow(qtldf))
      edf <- vector(length = nrow(qtldf))
      # determine x position of qtl
      
      if (nrow(qtldf) > 0) {
        for (ql in 1:nrow(qtldf)) {
          # if there are labels on y axis where qtl needs to be drawn
          # first decide if QTL line (so to eo) or text label is longest
          oneinch <- grconvertY(0, from = "user", to = "inches") -
            grconvertY(1, from = "user", to = "inches")
          lablen <- strwidth(qtldf$qtl[ql], units = "inches")
          lablenY <- lablen / oneinch
          if (lablenY > qtldf$eo[ql] - qtldf$so[ql]) {
            # use lablenY start is middle - half of label length
            qtlstart <-
              ((qtldf$si[ql] + qtldf$ei[ql]) / 2) - (lablenY / 2)
            qtlend <-
              ((qtldf$si[ql] + qtldf$ei[ql]) / 2) + (lablenY / 2)
          }
          else {
            qtlstart <- qtldf$so[ql]
            qtlend <- qtldf$eo[ql]
          }
          for (f in 1:length(sdf)) {
            if (qtlstart > sdf[f] && qtlstart < edf[f]) {
              qtlstart <- sdf[f]
            }
            if (qtlend < edf[f] && qtlend > sdf[f]) {
              qtlend <- edf[f]
            }
          }
          sdf[ql] <- qtlstart
          edf[ql] <- qtlend
          # adjust in case labels overlap
          
          
          if (any(adjyr >= (qtlstart + strheight("M") * max(rcex)) &
                  adjyr <= (qtlend   - strheight("M") * max(rcex)))) {
            qtlxpos <-
              max(yrlabwidth[yrlabwidth > 0][which(adjyr >= (qtlstart + strheight("M") * max(rcex)) &
                                                     adjyr <= (qtlend - strheight("M") * max(rcex)))]) / pdf.width
            yrlabwidth[yrlabwidth > 0][which(adjyr >= (qtlstart + strheight("M") * max(rcex)) &
                                               adjyr <= (qtlend  - strheight("M") * max(rcex)))] <-
              yrlabwidth[yrlabwidth > 0][which(adjyr >= (qtlstart + strheight("M") * max(rcex)) &
                                                 adjyr <= (qtlend  - strheight("M") * max(rcex)))] + lgw / 3 + strheight("M", units = "inches") *
              3
          }
          else {
            qtlxpos <- 0
          }
          
          if (is.null(qtldf$col[ql])) {
            qtldf$col[ql] <- par("col")
          }
          # draw qtl
          segments(
            x2[1] + posmult * (qtlxpos + strwidth("M") + lgwpct / 6),
            qtldf$so[ql],
            x2[1] + posmult * (qtlxpos + strwidth("M") + lgwpct / 6),
            qtldf$eo[ql],
            col = qtldf$col[ql]
          )
          # make inner region 1/3 size of linkage group width
          # 1/3 is arbitrary, just looks good
          rect(
            x2[1] + posmult * (qtlxpos + strwidth("M")),
            qtldf$si[ql],
            x2[1] + posmult * (qtlxpos + strwidth("M") + lgwpct / 3),
            qtldf$ei[ql],
            col = qtldf$col[ql],
            border = qtldf$col[ql]
          )
          text(
            x2[1] + posmult * (qtlxpos + strwidth("MM") + lgwpct / 3),
            (qtldf$si[ql] + qtldf$ei[ql]) / 2,
            label = qtldf$qtl[ql],
            # col = qtldf$col[ql],
            srt = 270
          )
        }
      }
    }
    return (
      list(
        adjyr = adjyr,
        adjyl = adjyr[dups$flkeep],
        yrlabwidth = returnlabwidth,
        dups = dups,
        solist = solist
      )
    )
  }

fsdups <- function(y, maxnbrcolsfordups) {
  k <- length(y)   # original length of y
  
  # For storing dup y locations
  yd <- matrix(nrow = maxnbrcolsfordups - 1, ncol = length(y))
  j <- 0           # dup counter
  
  lkeep <-
    vector(mode = "integer")  # keep track of which to keep labels for
  rkeep <-
    vector(mode = "integer")  # the rest go in columns no need to spread
  
  # keep a vector of indexes for the first column of duplicates
  # where more than one row of duplicates
  
  frkeep <- vector(mode = "integer")
  fykeep <- vector(mode = "integer")
  flkeep <- vector(mode = "integer")
  
  #always keep first one
  lkeep <- append(lkeep, 1)
  rkeep <- append(rkeep, 1)
  flkeep <- append(flkeep, 1)
  if (k > 1) {
    for (dupi in 2:k) {
      if (y[dupi] == y[dupi - 1]) {
        if (j == 0) {
          frkeep <- append(frkeep, dupi - 1)
          fykeep <- append(fykeep, length(rkeep))
        }
        if (j < maxnbrcolsfordups - 1) {
          j <- j + 1
          yd[j, (dupi - j)] <-
            length(rkeep)  # save what will be adjyr index for dup
        }
        else {
          j <- 0
          rkeep <- append(rkeep, dupi)
          frkeep <- append(frkeep, dupi)
          fykeep <-
            append(fykeep, length(rkeep))  # first one in subsequent row for dup
        }
      }
      else {
        if (j == 0 && dupi > 2 && y[(dupi - 1)] == y[(dupi - 2)]) {
          frkeep <- append(frkeep, dupi - 1)
          fykeep <- append(fykeep, length(rkeep))
        }
        j <- 0
        lkeep <- append(lkeep, dupi)
        rkeep <- append(rkeep, dupi)
        flkeep <-
          append(flkeep, length(rkeep))  # adjyr index for left labels
      }
    }
  }
  list(
    lkeep = lkeep[!duplicated(lkeep)],
    rkeep = rkeep[!duplicated(rkeep)],
    frkeep = frkeep[!duplicated(frkeep)],
    fykeep = fykeep[!duplicated(fykeep)],
    flkeep = flkeep[!duplicated(flkeep)],
    yd = yd
  )
}

lmvdencolor <- function(df,
                        wsize = 30,
                        bias = 5,
                        colorin = colorRampPalette(RColorBrewer::brewer.pal(8, "Spectral"))(25)) {
  editlgdf(df)
  mapthese <- unique(df[, 1])
  
  chr <- vector(mode = "character")
  s <- vector(mode = "double")
  e <- vector(mode = "double")
  col <- vector(mode = "character")
  dens <- vector(mode = "double")
  
  # assume if wsize=30 it defaulted.  If we need to adjust sectsize then we will
  # adjust wsize also since it defaulted.  If user passed it, use user value.
  if (wsize == 30) {
    adjwsize <- TRUE
  }
  else {
    adjwsize <- FALSE
  }
  
  # determine size of section to color.  Increase size until < 1000 sections for
  # largest linkage group (in case input is in Mbp instead of cM)
  sectsize <- 1
  if (max(df[, 2]) > 1000) {
    nbrsects <- max(df[, 2])
    sectsize <- 1
    while (nbrsects > 1000) {
      sectsize <- sectsize * 10
      # adjust wsize in the same way
      if (adjwsize) {
        wsize <- wsize * 10
      }
      nbrsects <- max(df[, 2]) / sectsize
    }
  }
  # calculate density within linkage group in sectsize position units
  dfix <- 0
  for (i in 1:length(mapthese)) {
    thislg <- df[, 2][df[, 1] == mapthese[i]]
    k <- min(thislg) + sectsize / 2
    while (k <= max(thislg) - sectsize / 2) {
      dfix <- dfix + 1
      chr[dfix] <- mapthese[i]
      s[dfix] <- k - sectsize / 2
      # if < 1 sectionleft to end go to end
      if ((k + (sectsize/2 * 3)) > max(thislg)) {
        e[dfix] <- max(thislg)
      }
      else {
        e[dfix] <- k + sectsize / 2
      }
      if (sum(thislg < (k + wsize / 2) &
              thislg > (k - wsize / 2)) > 0) {
        # calculate actual window size
        awsize <-
          min(k + wsize / 2, max(thislg)) - max(k - wsize / 2, min(thislg))
        dens[dfix] <- awsize / (sum(thislg <= (k + wsize / 2) &
                                      thislg >= (k - wsize / 2)))
      }
      else {
        # means we have a window size region with no markers
        # so cM/locus is infinity
        message(c("There are sections of linkage group ", mapthese[i],
                  " with no markers found between ",
                  k - wsize / 2, " and ",  k + wsize / 2,". Your window size is ", wsize, "."))
        message("These sections will be colored at the highest density.")
        message("Or use function lmvdencolor to set an alternative window size. ")
        dens[dfix] <- NA  # will set to 1 after scaling
      }
      col[dfix] = "default"
      k <- k + sectsize
    }
  }
  sectcoldf <-
    data.frame(chr, s, e, col, dens, stringsAsFactors = FALSE)
  
  # scale all densities to range of 0 to 1
  densrange <- range(sectcoldf$dens, na.rm = TRUE) #ignore NAs
  if (diff(densrange) > 0) {
    densvals  <- (sectcoldf$dens - densrange[1]) / diff(densrange)
  }
  else
  {
    densvals <- rep(0, length(sectcoldf$dens))
  }
  # set NA to 1 as highest position / marker
  densvals[is.na(densvals)] <- 1.0
  
  f <- colorRamp(colorin, bias = bias)  # create color ramp function
  
  sectcoldf$col <- rgb(f(densvals) / 255)
  sectcoldf
  
}

readlgcross <- function(cross, mapthese)
{
  
  if (!requireNamespace("qtl", quietly = TRUE)) {
    stop("Please install package qtl")
  }
  # if list not provided create list of all using chrnames function of qtl
  if (missing(mapthese)) {
    mapthese <- qtl::chrnames(cross)
  }
  else {
    if (!all(mapthese %in% qtl::chrnames(cross))) {
      stop ("chrnames to map not found in cross object")
    }
  }
  df <- qtl::pull.map(cross, mapthese, as.table = TRUE)
  
  # remove factors if any
  fas <- sapply(df, is.factor)
  df[fas] <- lapply(df[fas], as.character)
  
  df$locus <- rownames(df)
  
  # assign my own column names so I can reference by name
  colnames(df) <- c("group", "position", "locus")
  return (df)
}

readlgtext <-
  function(fn,
           mapthese,
           header = TRUE,
           stringsAsFactors = FALSE)
  {
    # get file extension
    fnparts <- strsplit(fn, ".", fixed = TRUE)[[1]]
    if (length(fnparts) > 1) {
      ext = fnparts[length(fnparts)]
    } else {
      ext <- ""
    }
    
    # figure out how many rows and columns on input file
    # and issue message
    
    if (ext == 'csv') {
      nbrcols <- count.fields(fn, sep = ",")
    }
    else {
      nbrcols <- count.fields(fn)
    }
    
    if (any(nbrcols < 2)) {
      for (i in 1:length(nbrcols)) {
        if (nbrcols[i] < 2) {
          message(c("line ", i, " has only ", nbrcols[i], " columns."))
        }
      }
      stop("Less than 2 columns found on input file")
    }
    if (any(nbrcols < 3)) {
      for (i in 1:length(nbrcols)) {
        if (nbrcols[i] < 3) {
          message(c("row ", i, " is missing locus name."))
        }
      }
      warning(" Locus name missing on some lines")
    }
    
    if (ext == 'csv') {
      df <- read.csv(
        fn,
        header = header,
        stringsAsFactors = FALSE,
        strip.white = TRUE
      )
    }
    else {
      df <-
        read.table(
          fn,
          header = header,
          stringsAsFactors = FALSE,
          colClasses = c("character", "numeric", "character"),
          strip.white = TRUE
        )
      
    }
    # assign my own column names so I can reference by name
    colnames(df)[1:3] <- c("group", "position", "locus")
    
    # make group a character field
    
    df$group <- as.character(df$group)
    
    if (!missing(mapthese)) {
      if (!all(mapthese %in% df$group)) {
        stop ("chrnames to map not found in input file")
      }
    }
    
    return (df)
  }

revpos <- function(pos, maxdec)
{
  newpos <- vector(length = length(pos))
  newpos[1] <- pos[1]
  for (i in 2:length(pos)) {
    newpos[i] <-
      newpos[(i - 1)] + pos[(length(pos) - i + 2)] - pos[(length(pos) - i + 1)]
  }
  return (round(newpos, digits = maxdec))
}

show <- function(showonly,
                 llab,
                 rlab,
                 rcex = par("cex"),
                 lcex = par("cex"),
                 rfont = par("font"),
                 lfont = par("font"),
                 rcol = par("col"),
                 lcol = par("col"))
{
  # when invoked all vectors except showonly should be the
  # same length - see lmv.linkgae.plot.R where they are set
  
  newllab <- vector()
  newrlab <- vector()
  newrcex <- vector()
  newlcex <- vector()
  newrfont <- vector()
  newlfont <- vector()
  newrcol <- vector()
  newlcol <- vector()
  
  # match(showonly, rlab) returns the index where there is a match
  # since the result is in showonly order, the index in rlab might
  # not be in order, sort them.  that also removes the NA since the
  # showonly might not be in this linkage group
  
  newllab <- append(newllab, llab[sort(match(showonly, rlab))])
  newrlab <- append(newrlab, rlab[sort(match(showonly, rlab))])
  newrcex <- append(newrcex, rcex[sort(match(showonly, rlab))])
  newlcex <- append(newlcex, lcex[sort(match(showonly, rlab))])
  newrfont <- append(newrfont, rfont[sort(match(showonly, rlab))])
  newlfont <- append(newlfont, lfont[sort(match(showonly, rlab))])
  newrcol <- append(newrcol, rcol[sort(match(showonly, rlab))])
  newlcol <- append(newlcol, lcol[sort(match(showonly, rlab))])
  
  
  # always show the first and last position label
  
  if (!(rlab[1] %in% showonly)) {
    newrlab <- append("", newrlab)
    newllab <- append(llab[1], newllab)
    newrcex <- append(rcex[1], newrcex)
    newlcex <- append(lcex[1], newlcex)
    newrfont <- append(rfont[1], newrfont)
    newlfont <- append(lfont[1], newlfont)
    newrcol <- append(rcol[1], newrcol)
    newlcol <- append(lcol[1], newlcol)
  }
  
  if (!(rlab[length(rlab)] %in% showonly)) {
    newrlab <- append(newrlab, "")
    newllab <- append(newllab, llab[length(llab)])
    newrcex <- append(newrcex, rcex[length(llab)])
    newlcex <- append(newlcex, lcex[length(llab)])
    newrfont <- append(newrfont, rfont[length(llab)])
    newlfont <- append(newlfont, lfont[length(llab)])
    newrcol <- append(newrcol, rcol[length(llab)])
    newlcol <- append(newlcol, lcol[length(llab)])
  }
  
  return (
    list(
      newllab = newllab,
      newrlab = newrlab,
      newrcex = newrcex,
      newlcex = newlcex,
      newrfont = newrfont,
      newlfont = newlfont,
      newrcol = newrcol,
      newlcol = newlcol
    )
  )
}

spreadcexlabs <-
  function(x,
           mindiff,
           strh,
           maxiter = 7500,
           stepsize = 1 / 10,
           min = -Inf,
           max = Inf,
           cex = par("cex")) {
    # assumes ordered input due to variable sizes of labels based on cex
    
    cex = rep(cex, length.out = length(x))
    mingap <- mindiff - strh
    
    reqdiff <- vector(length = length(x) - 1)
    
    reqdiff[1:(length(x) - 1)] <-
      0.5 * strh * cex[1:(length(x) - 1)] + 0.5 * strh * cex[2:(length(x))] + mingap
    
    df <-
      x[-1] - x[-length(x)]  # get differences between points
    
    stp <- mindiff * stepsize
    
    i <- 1
    while (any(df < reqdiff)) {
      tmp <- c(df < reqdiff, FALSE)
      if (tmp[1] &&
          (x[1] - stp) < min) {
        # don't move bottom set
        tmp2 <- as.logical(cumprod(tmp))
        tmp <- tmp & !tmp2
      }
      x[tmp] <- x[tmp] - stp  # move all non-fitters down
      
      
      tmp <- c(FALSE, df < reqdiff)
      if (tmp[length(tmp)] &&
          (x[length(x)] + stp) > max) {
        # don't move top
        tmp2 <- rev(as.logical(cumprod(rev(tmp))))
        tmp <- tmp & !tmp2
      }
      x[tmp] <- x[tmp] + stp   # move new all non-fitters up
      
      df <- x[-1] - x[-length(x)]
      i <- i + 1
      if (i > maxiter) {
        #                warning("Maximum iterations reached")
        break
      }
    }
    
    return(x)
  }

usescanone <- function(qtlscanone, qtldf=NULL, mapthese, fg, maxdec) {
  
  if (is.null(qtldf)) {
    # make a dataframe to hold qtl info
    qtldf <- data.frame(
      chr = character(),
      qtl = character(),
      so = numeric(),
      si = numeric(),
      ei = numeric(),
      eo = numeric(),
      col = character(),
      stringsAsFactors = FALSE
    )
  }
  # extract qtl info from scanone df
  qtlchr <- unique(qtlscanone[, 1])
  for (i in 1:length(qtlchr)) {
    if (qtlchr[i] %in% mapthese) {
      qtlone <- qtl::bayesint(qtlscanone, chr = qtlchr[i])
      qtldf <- rbind(
        qtldf,
        data.frame(
          chr = qtlchr[i],
          qtl = paste(qtlchr[i], "@" , round(qtlone[2,2], digits = maxdec)),
          so = qtlone[1, 2],
          si = qtlone[2, 2],
          ei = qtlone[2, 2],
          eo = qtlone[3, 2],
          col=fg
        )
      )
    }
  }
  # remove factors if any
  fas <- sapply(qtldf, is.factor)
  qtldf[fas] <- lapply(qtldf[fas], as.character)
  return (qtldf)
}