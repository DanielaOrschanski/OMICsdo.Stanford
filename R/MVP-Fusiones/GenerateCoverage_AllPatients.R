library(GenomicRanges)
library(GenomicAlignments)


patients_dir <- "/media/16TBDisk/Daniela/Fusiones/Gastrico-610"
ids <- list.dirs(patients_dir, recursive = FALSE, full.names = FALSE)

# DataFrame de fusiones detectadas
# Ejemplo simplificado: deberías tener uno por paciente
# con columnas gene_id1, contig1, breakpoint1, gene_id2, contig2, breakpoint2
ruta_archivo <- "/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/Act_ReportFusions_TodasCohortes.xlsx"
fusions <- as.data.frame(read_excel(ruta_archivo))
#fusions <- as.data.frame(fusions_df[which(fusions_df$ID %in% ids),])
fusions <- fusions[which(fusions$ID == "39738"),]
unique(fusions$ID)
fusions <- fusions[2,]

load("/media/16TBDisk/Daniela/OMICsdo/Grafico-Fusiones-Propio/exons-arriba.RData")
source("/media/16TBDisk/Daniela/OMICsdo/Grafico-Fusiones-Propio/draw_fusions_PROPIO.R", echo=TRUE)
plot_Fusion_without_coverage(fusions, exons)

plot_Fusion_without_coverage <- function(fusions, exons) {

  fusion = 1
  sampleName <- fusions$ID
  message(paste0("Drawing fusion #", fusion, ": ", fusions[fusion,"gene1"], ":", fusions[fusion,"gene2"]))

  if(fusions$peptide_sequence[fusion] == ".") {
    stop( message(paste0("Fusion: ", fusions[fusion,"gene1"], ":", fusions[fusion,"gene2"], "did not generate a protein.")))
  }

  # if showIntergenicVicinity is a number, take it as is
  # if it is a keyword (closestGene/closestProteinCodingGene), determine the range dynamically
  showVicinity <- rep(0, 4)
  if (fusions[fusion,"site1"] == "intergenic") {
    showVicinity[1] <- ifelse(
      is.numeric(showIntergenicVicinity[[1]]),
      showIntergenicVicinity[[1]],
      fusions[fusion,"breakpoint1"] - start(findClosestGene(exons, fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"], exons$end < fusions[fusion,"breakpoint1"] & exons$type == showIntergenicVicinity[[1]]))
    )
    showVicinity[2] <- ifelse(
      is.numeric(showIntergenicVicinity[[2]]),
      showIntergenicVicinity[[2]],
      end(findClosestGene(exons, fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"], exons$start > fusions[fusion,"breakpoint1"] & exons$type == showIntergenicVicinity[[2]])) - fusions[fusion,"breakpoint1"]
    )
  }
  if (fusions[fusion,"site2"] == "intergenic") {
    showVicinity[3] <- ifelse(
      is.numeric(showIntergenicVicinity[[3]]),
      showIntergenicVicinity[[3]],
      fusions[fusion,"breakpoint2"] - start(findClosestGene(exons, fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"], exons$end < fusions[fusion,"breakpoint2"] & exons$type == showIntergenicVicinity[[3]]))
    )
    showVicinity[4] <- ifelse(
      is.numeric(showIntergenicVicinity[[4]]),
      showIntergenicVicinity[[4]],
      end(findClosestGene(exons, fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"], exons$start > fusions[fusion,"breakpoint2"] & exons$type == showIntergenicVicinity[[4]])) - fusions[fusion,"breakpoint2"]
    )
  }

  # compute coverage from alignments file
  # find all exons belonging to the fused genes
  coverage1 <- NULL
  coverage2 <- NULL
  exons1 <- findExons(exons, fusions[fusion,"contig1"], fusions[fusion,"gene_id1"], fusions[fusion,"direction1"], fusions[fusion,"breakpoint1"], coverage1, fusions[fusion,"transcript_id1"], transcriptSelection)
  if (nrow(exons1) == 0) {
    par(mfrow=c(1,1))
    plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
    text(0, 0, paste0("exon coordinates of ", fusions[fusion,"gene1"], " not found in\n", exonsFile))
    warning(paste("exon coordinates of", fusions[fusion,"gene1"], "not found"))
    next
  }
  exons2 <- findExons(exons, fusions[fusion,"contig2"], fusions[fusion,"gene_id2"], fusions[fusion,"direction2"], fusions[fusion,"breakpoint2"], coverage2, fusions[fusion,"transcript_id2"], transcriptSelection)
  if (nrow(exons2) == 0) {
    par(mfrow=c(1,1))
    plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
    text(0, 0, paste0("exon coordinates of ", fusions[fusion,"gene2"], " not found in\n", exonsFile))
    warning(paste("exon coordinates of", fusions[fusion,"gene2"], "not found"))
    next
  }

  # in case of intergenic breakpoints, show the vicinity
  if (sum(showVicinity) > 0) {
    if (fusions[fusion,"site1"] == "intergenic") {
      for (geneID in unique(exons[exons$contig == fusions[fusion,"contig1"] & exons$exonNumber != "intergenic" &
                                  (between(exons$end, fusions[fusion,"breakpoint1"]-showVicinity[1], fusions[fusion,"breakpoint1"]+showVicinity[2]) |
                                   between(exons$start, fusions[fusion,"breakpoint1"]-showVicinity[1], fusions[fusion,"breakpoint1"]+showVicinity[2])),"geneID"]))
        exons1 <- rbind(exons1, findExons(exons, fusions[fusion,"contig1"], geneID, fusions[fusion,"direction1"], fusions[fusion,"breakpoint1"], coverage1, fusions[fusion,"transcript_id1"], transcriptSelection))
      # crop genes that are only partially within user-defined vicinity, because the coverage data is incomplete for those
      exons1 <- exons1[exons1$start >= fusions[fusion,"breakpoint1"]-showVicinity[1] & exons1$end <= fusions[fusion,"breakpoint1"]+showVicinity[2] | exons1$exonNumber == "intergenic",]
    }
    if (fusions[fusion,"site2"] == "intergenic") {
      for (geneID in unique(exons[exons$contig == fusions[fusion,"contig2"] & exons$exonNumber != "intergenic" &
                                  (between(exons$end, fusions[fusion,"breakpoint2"]-showVicinity[3], fusions[fusion,"breakpoint2"]+showVicinity[4]) |
                                   between(exons$start, fusions[fusion,"breakpoint2"]-showVicinity[3], fusions[fusion,"breakpoint2"]+showVicinity[4])),"geneID"]))
        exons2 <- rbind(exons2, findExons(exons, fusions[fusion,"contig2"], geneID, fusions[fusion,"direction2"], fusions[fusion,"breakpoint2"], coverage2, fusions[fusion,"transcript_id2"], transcriptSelection))
      # crop genes that are only partially within user-defined vicinity, because the coverage data is incomplete for those
      exons2 <- exons2[exons2$start >= fusions[fusion,"breakpoint2"]-showVicinity[3] & exons2$end <= fusions[fusion,"breakpoint2"]+showVicinity[4] | exons2$exonNumber == "intergenic",]
    }
  }

  # sort coding exons last, such that they are drawn over the border of non-coding exons
  exons1 <- exons1[order(exons1$start, -rank(exons1$type)),]
  exons2 <- exons2[order(exons2$start, -rank(exons2$type)),]

  # insert dummy exons, if breakpoints are outside the gene (e.g., in UTRs)
  # this avoids plotting artifacts
  breakpoint1 <- fusions[fusion,"breakpoint1"]
  breakpoint2 <- fusions[fusion,"breakpoint2"]
  if (breakpoint1 < min(exons1$start)) {
    exons1 <- rbind(c(exons1[1,"contig"], "dummy", max(1,breakpoint1-1000), max(1,breakpoint1-1000), exons1[1,"strand"], "", "dummy", exons1[1,"geneID"], exons1[1,"transcript"], ""), exons1)
  } else if (breakpoint1 > max(exons1$end)) {
    exons1 <- rbind(exons1, c(exons1[1,"contig"], "dummy", breakpoint1+1000, breakpoint1+1000, exons1[1,"strand"], "", "dummy", exons1[1,"geneID"], exons1[1,"transcript"], ""))
  }
  if (breakpoint2 < min(exons2$start)) {
    exons2 <- rbind(c(exons2[1,"contig"], "dummy", max(1,breakpoint2-1000), max(1,breakpoint2-1000), exons2[1,"strand"], "", "dummy", exons2[1,"geneID"], exons2[1,"transcript"], ""), exons2)
  } else if (breakpoint2 > max(exons2$end)) {
    exons2 <- rbind(exons2, c(exons2[1,"contig"], "dummy", breakpoint2+1000, breakpoint2+1000, exons2[1,"strand"], "", "dummy", exons2[1,"geneID"], exons2[1,"transcript"], ""))
  }
  exons1$start <- as.integer(exons1$start)
  exons1$end <- as.integer(exons1$end)
  exons2$start <- as.integer(exons2$start)
  exons2$end <- as.integer(exons2$end)

  exons1$left <- exons1$start
  exons1$right <- exons1$end
  exons2$left <- exons2$start
  exons2$right <- exons2$end

  squishedIntronSize <- 200
  if (squishIntrons) {
    # hide introns in gene1
    cumulativeIntronLength <- 0
    previousExonEnd <- -squishedIntronSize
    for (exon in 1:nrow(exons1)) {
      if (breakpoint1 > previousExonEnd+1 && breakpoint1 < exons1[exon,"left"])
        breakpoint1 <- (breakpoint1-previousExonEnd) / (exons1[exon,"left"]-previousExonEnd) * squishedIntronSize + previousExonEnd - cumulativeIntronLength
      if (exons1[exon,"left"] > previousExonEnd) {
        cumulativeIntronLength <- cumulativeIntronLength + exons1[exon,"left"] - previousExonEnd - squishedIntronSize
        previousExonEnd <- exons1[exon,"right"]
      }
      if (breakpoint1 >= exons1[exon,"left"] && breakpoint1 <= exons1[exon,"right"]+1)
        breakpoint1 <- breakpoint1 - cumulativeIntronLength
      exons1[exon,"left"] <- exons1[exon,"left"] - cumulativeIntronLength
      exons1[exon,"right"] <- exons1[exon,"right"] - cumulativeIntronLength
    }

    # hide introns in gene2
    cumulativeIntronLength <- 0
    previousExonEnd <- -squishedIntronSize
    for (exon in 1:nrow(exons2)) {
      if (breakpoint2 > previousExonEnd+1 && breakpoint2 < exons2[exon,"left"])
        breakpoint2 <- (breakpoint2-previousExonEnd) / (exons2[exon,"left"]-previousExonEnd) * squishedIntronSize + previousExonEnd - cumulativeIntronLength
      if (exons2[exon,"left"] > previousExonEnd) {
        cumulativeIntronLength <- cumulativeIntronLength + exons2[exon,"left"] - previousExonEnd - squishedIntronSize
        previousExonEnd <- exons2[exon,"right"]
      }
      if (breakpoint2 >= exons2[exon,"left"] && breakpoint2 <= exons2[exon,"right"]+1)
        breakpoint2 <- breakpoint2 - cumulativeIntronLength
      exons2[exon,"left"] <- exons2[exon,"left"] - cumulativeIntronLength
      exons2[exon,"right"] <- exons2[exon,"right"] - cumulativeIntronLength
    }
  } else { # don't squish introns
    # shift exon coordinates to align the gene to the left border of the plot
    exons1$right <- exons1$right - min(exons1$left)
    breakpoint1 <- breakpoint1 - min(exons1$left)
    exons1$left <- exons1$left - min(exons1$left)
    exons2$right <- exons2$right - min(exons2$left)
    breakpoint2 <- breakpoint2 - min(exons2$left)
    exons2$left <- exons2$left - min(exons2$left)
  }

  # scale exon sizes to fit on page
  scalingFactor <- max(exons1$right) + max(exons2$right)
  if (fixedScale > 0) {
    if (fixedScale >= scalingFactor) {
      scalingFactor <- fixedScale
    } else {
      warning(paste("fallback to automatic scaling, because value for --fixedScale is too small to fit transcripts on canvas (increase it to", scalingFactor, "to avoid this)"))
    }
  }
  exons1$left <- exons1$left / scalingFactor
  exons1$right <- exons1$right / scalingFactor
  exons2$left <- exons2$left / scalingFactor
  exons2$right <- exons2$right / scalingFactor
  breakpoint1 <- breakpoint1 / scalingFactor
  breakpoint2 <- breakpoint2 / scalingFactor

  # shift gene2 to the right border of the page
  #gene2Offset <- 1 + 0.05 - max(exons2$right)

  gene2Offset <- 0.2 - max(exons2$right)

  # center fusion horizontally
  fusionOffset1 <- (max(exons1$right)+gene2Offset)/2 - ifelse(fusions[fusion,"direction1"] == "downstream", breakpoint1, max(exons1$right)-breakpoint1)
  fusionOffset1 <- fusionOffset1 + 0.4
  fusionOffset2 <- fusionOffset1 + ifelse(fusions[fusion,"direction1"] == "downstream", breakpoint1, max(exons1$right)-breakpoint1)
  #fusionOffset2 <- fusionOffset2 + 0.2

  ########################################
  # ARRANCA GRAFICO ########################################
  ####################################

  # layout: fusion on top, circos plot on bottom left, protein domains on bottom center, statistics on bottom right
  layout(matrix(c(1,1,1,2,4,5,3,4,5), 3, 3, byrow=TRUE), widths=c(1.1, 1.2, 0.7), heights=c(1.55, 1.2, 0.25))
  par(mar=c(0, 0, 0, 0))
  plot(0, 0, type="l", xlim=c(-0.12, 1.12), ylim=c(0, 1.1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  fontSize = 1.5
  # vertical coordinates of layers
  ySampleName <- 1
  yBreakpointLabels <- 0.7
  yGeneNames <- 0.9
  yFusion <- 0.6
  yTranscript <- 0.5
  yScale <- 0.407

  # print sample name (title of page)
  text(0.5, ySampleName, sampleName, font=2, cex=fontSize*1.2, adj=c(0.5,0))


  # draw gene & transcript names
  #fusion = 2
  if (fusions[fusion,"gene1"] != ".")
    #text(max(exons1$right)/2, yGeneNames, fusions[fusion,"gene1"], font=2, cex=fontSize, adj=c(0.5,0))
    text(0.25, yGeneNames, fusions[fusion,"gene1"], font=2, cex=fontSize, adj=c(0.5,0))
  if (fusions[fusion,"site1"] != "intergenic")
    #text(max(exons1$right)/2, yGeneNames-0.01, head(exons1$transcript,1), cex=0.9*fontSize, adj=c(0.5,1))
    text(0.25, yGeneNames-0.02, sprintf("%s", head(exons1$transcript,1)), cex=0.9*fontSize, adj=c(0.5,1))

  if (fusions[fusion,"gene2"] != ".")
    #text(gene2Offset+max(exons2$right)/2, yGeneNames, fusions[fusion,"gene2"], font=2, cex=fontSize, adj=c(0.5,0))
    text(0.75, yGeneNames, fusions[fusion,"gene2"], font=2, cex=fontSize, adj=c(0.5,0))
  if (fusions[fusion,"site2"] != "intergenic")
    #text(gene2Offset+max(exons2$right)/2, yGeneNames-0.01, head(exons2$transcript,1), cex=0.9*fontSize, adj=c(0.5,1))
    text(0.75, yGeneNames-0.02, sprintf("%s", head(exons2$transcript,1)), cex=0.9*fontSize, adj=c(0.5,1))

  # if multiple genes in the vicinity are shown, label them
  if (fusions[fusion,"site1"] == "intergenic")
    for (gene in unique(exons1$geneName)) {
      exonsOfGene <- exons1[exons1$geneName == gene & exons1$type != "dummy",]
      if (any(exonsOfGene$type == "exon"))
        text(mean(c(min(exonsOfGene$left), max(exonsOfGene$right))), yGeneNames-0.04, gene, cex=0.9*fontSize, adj=c(0.5,1))
    }
  if (fusions[fusion,"site2"] == "intergenic")
    for (gene in unique(exons2$geneName)) {
      exonsOfGene <- exons2[exons2$geneName == gene & exons2$type != "dummy",]
      if (any(exonsOfGene$type == "exon"))
        text(gene2Offset+mean(c(min(exonsOfGene$left), max(exonsOfGene$right))), yGeneNames-0.04, gene, cex=0.9*fontSize, adj=c(0.5,1))
    }

  # label breakpoints
  yBreakpointLabels = 0.75
  #text(breakpoint1+0.01, yBreakpointLabels-0.03, paste0("bp1 = ", fusions[fusion,"display_contig1"], ":", fusions[fusion,"breakpoint1"]), adj=c(1,0), cex=fontSize)
  text(0.25, yBreakpointLabels, paste0("bp1 = chr ", fusions[fusion,"display_contig1"], ":", fusions[fusion,"breakpoint1"]), adj=c(0.5,0), cex=fontSize*0.9)
  #text(gene2Offset+breakpoint2-0.01, yBreakpointLabels-0.03, paste0("breakpoint2\n", fusions[fusion,"display_contig2"], ":", fusions[fusion,"breakpoint2"]), adj=c(0,0), cex=fontSize)
  text(0.75, yBreakpointLabels, paste0("bp2 = chr ", fusions[fusion,"display_contig2"], ":", fusions[fusion,"breakpoint2"]), adj=c(0.5,0), cex=fontSize*0.9)

  # draw coverage axis
  exons1 <- exons1[exons1$type == "exon",]
  exons2 <- exons2[exons2$type == "exon",]

  #-----------------------------------------------------------------------

  # plot gene1 of fusion -----------
  margin <- 0.05

  if (fusions[fusion,"direction1"] == "downstream") {
    # plot strands
    lines(c(fusionOffset1, fusionOffset1+breakpoint1), c(yFusion, yFusion), col=darkColor1)
    for (gene in unique(exons1$geneName)) {
      exonsOfGene <- exons1[exons1$geneName == gene,]
      if (min(exonsOfGene$start) <= fusions[fusion,"breakpoint1"])
        drawStrand(fusionOffset1+min(exonsOfGene$left), fusionOffset1+min(breakpoint1, max(exonsOfGene$right)), yFusion, col=darkColor1, exonsOfGene$strand[1])
    }
    # plot exons
    for (exon in 1:nrow(exons1))
      if (exons1[exon,"start"] <= fusions[fusion,"breakpoint1"])
        drawExon_bigger(fusionOffset1+exons1[exon,"left"], fusionOffset1+min(breakpoint1, exons1[exon,"right"]), yFusion, color = color1, color_claro = "#DCD0FF", exons1[exon,"exonNumber"], exons1[exon,"type"])

  } else if (fusions[fusion,"direction1"] == "upstream") {

    for (gene in unique(exons1$geneName)) {
      exonsOfGene <- exons1[exons1$geneName == gene,]
      if (max(exonsOfGene$end+1) >= fusions[fusion,"breakpoint1"])
        drawStrand(breakpoint1, max(exonsOfGene$right), yFusion, col=darkColor1, chartr("+-", "-+", exonsOfGene$strand[1]))
    }
    # plot exons
    for (exon in 1:nrow(exons1))
      if (exons1[exon,"end"]+1 >= fusions[fusion,"breakpoint1"])
        drawExon_bigger(fusionOffset1+max(exons1$right)-exons1[exon,"right"], min(fusionOffset2, fusionOffset1+max(exons1$right)-exons1[exon,"left"]), yFusion, color = color1, color_claro = "#DCD0FF", exons1[exon,"exonNumber"], exons1[exon,"type"])
  }


  # plot gene2 of fusion ---
  if (fusions[fusion,"direction2"] == "downstream") {
    # plot strands
    lines(c(fusionOffset2, fusionOffset2+breakpoint2), c(yFusion, yFusion), col=darkColor2)
    for (gene in unique(exons2$geneName)) {
      exonsOfGene <- exons2[exons2$geneName == gene,]
      if (min(exonsOfGene$start) <= fusions[fusion,"breakpoint2"])
        drawStrand(max(fusionOffset2, fusionOffset2+breakpoint2-max(exonsOfGene$right)), fusionOffset2+breakpoint2-min(exonsOfGene$left), yFusion, col=darkColor2, chartr("+-", "-+", exonsOfGene$strand[1]))
    }
    # plot exons
    for (exon in 1:nrow(exons2))
      if (exons2[exon,"start"] <= fusions[fusion,"breakpoint2"])
        drawExon_bigger(max(fusionOffset2, fusionOffset2+breakpoint2-exons2[exon,"right"]), fusionOffset2+breakpoint2-exons2[exon,"left"], yFusion, color = color2, color_claro = "#FFDBBB", exons2[exon,"exonNumber"], exons2[exon,"type"])
  } else if (fusions[fusion,"direction2"] == "upstream") {
    # plot strands
    lines(c(fusionOffset2, fusionOffset2+max(exons2$right)-breakpoint2), c(yFusion, yFusion), col=darkColor2)
    for (gene in unique(exons2$geneName)) {
      exonsOfGene <- exons2[exons2$geneName == gene,]
      if (max(exonsOfGene$end+1) >= fusions[fusion,"breakpoint2"])
        drawStrand(max(fusionOffset2, fusionOffset2+min(exonsOfGene$left)-breakpoint2), fusionOffset2+max(exonsOfGene$right)-breakpoint2, yFusion, col=darkColor2, exonsOfGene$strand[1])
    }
    # plot exons
    for (exon in 1:nrow(exons2))
      if (exons2[exon,"end"]+1 >= fusions[fusion,"breakpoint2"])
        drawExon_bigger(max(fusionOffset2, fusionOffset2+exons2[exon,"left"]-breakpoint2), fusionOffset2+exons2[exon,"right"]-breakpoint2, yFusion, color = color2, color_claro = "#FFDBBB", exons2[exon,"exonNumber"], exons2[exon,"type"])
  }
  # print statistics about supporting alignments: ---------------------------
  text(fusionOffset2-0.1, 0.35, "Supporting Read Count", font=2, adj=c(0,0), cex=1)
  text(
    fusionOffset2-0.1, 0.325,
    paste0(
      "Split reads - bp1 = ", fusions[fusion,"split_reads1"], "\n",
      "Split reads - bp2 = ", fusions[fusion,"split_reads2"], "\n",
      "Discordant mates = ", fusions[fusion,"discordant_mates"], "\n",
      "Confidence = ", fusions[fusion, "confidence"]
    ),
    adj=c(0,1), cex=1
  )
}

