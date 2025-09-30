#' @title downloadPfam
#' @description Generates
#' @param omicsdo_sof_path path
#' @export
#' @import openxlsx
#' @import readr
downloadPfam <- function(omicsdo_sof_path) {
  setwd(omicsdo_sof_path)

  pfam_path <- sprintf("%s/Pfam-A.hmm", omicsdo_sof_path)

  if(file.exists(pfam_path)) {
    message("Pfam was already downloaded")
  } else {
    setwd(dirname(pfam_path))
    message("Pfam database will be downloaded")
    system("wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.hmm.gz")
    system(sprintf("gunzip %s/Pfam-A.hmm.gz", dirname(pfam_path)))
    #system2("sudo apt install hmmer") # hacer esto en el terminal
    system(sprintf("hmmpress %s ", pfam_path))
  }
  return(pfam_path)
}

# Execution:
#protein_df <- readRDS("/media/respaldo8t/Daniela/Transcript_PeptideSequence.RDS")
#fusions_report <- read_excel("/media/8tb02/Daniela/Leucemia/Todos-FusionReports_Leucemia.xlsx")
#colnames(fusions_report)[1] <- "ID"

#fusions_report_anotated <- AnotateFusions(protein_df, fusions_report)
#fusions_report_anotated_domain <- GetDomains(fusions_report = fusions_report_anotated, path_dir = path_dir)
#domain = domain_1
#fusion = seq_prot


get_longest_common_substring <- function(domain, fusion) {
  domain_seq <- AAString(domain)
  fusion_seq <- AAString(fusion)

  hits <- pairwiseAlignment(domain_seq, fusion_seq, type = "local", scoreOnly = FALSE)
  conserved_domain <- as.character(subject(hits))
  return(conserved_domain)
}

#' @title GetDomains
#' @description Generates
#' @param fusions_report df
#' @param path_dir path to
#' @export
#' @import openxlsx
#' @import readr


#domains <- GetDomains(gene1 = "KIAA1549",
#                      gene2 = "BRAF",
#                      transcript_id1 = "ENST00000422774",
#                      transcript_id2 = "ENST00000496384",
#                      seq_peptide1 = "MPGARRRRRGAAMEGKPRAGVALAPGPSGRRPSARCARRRRPGLLLPGLWLLLLARPASCAPDELSPEQHNLSLYSMELVLKKSTGHSAAQVALTETAPGSQHSSPLHVTAPPSATTFDTAFFNQGKQTKSTADPSIFVATYVSVTSKEVAVNDDEMDNFLPDTHWTTPRMVSPIQYITVSPPGLPREALEPMLTPSLPMVSLQDEEVTSGWQNTTRQPAAYAESASHFHTFRSAFRTSEGIVPTPGRNLVLYPTDAYSHLSSRTLPEIVASLTEGVETTLFLSSRSLMPQPLGDGITIPLPSLGEVSQPPEEVWATSADRYTDVTTVLSQSLEETISPRTYPTVTASHAALAFSRTHSPLLSTPLAFASSASPTDVSSNPFLPSDSSKTSELHSNSALPGPVDNTHILSPVSSFRPYTWCAACTVPSPQQVLATSLMEKDVGSGDGAETLCMTVLEESSISLMSSVVADFSEFEEDPQVFNTLFPSRPIVPLSSRSMEISETSVGISAEVDMSSVTTTQVPPAHGRLSVPASLDPTAGSLSVAETQVTPSSVTTAFFSVITSILLDSSFSVIANKNTPSLAVRDPSVFTPYSLVPSVESSLFSDQERSSFSEHKPRGALDFASSFFSTPPLELSGSISSPSEAPASLSLMPSDLSPFTSQSFSPLVETFTLFDSSDLQSSQLSLPSSTNLEFSQLQPSSELPLNTIMLLPSRSEVSPWSSFPSDSLEFVEASTVSLTDSEAHFTSAFIETTSYLESSLISHESAVTALVPPGSESFDILTAGIQATSPLTTVHTTPILTESSLFSTLTPPDDQISALDGHVSVLASFSKAIPTGTVLITDAYLPSGSSFVSEATPFPLPTELTVVGPSLTPTEVPLNTSTEVSTTSTGAATGGPLDSTLMGDAASQSPPESSAAPPLPSLRPVTAFTLEATVDTPTLATAKPPYVCDITVPDAYLITTVLARRAVQEYIITAIKEVLRIHFNRAVELKVYELFTDFTFLVTSGPFVYTAISVINVLINSKLVRDQTPLILSVKPSFLVPESRFQVQTVLQFVPPSVDTGFCNFTQRIEKGLMTALFEVRKHHQGTYNLTVQILNITISSSRVTPRRGPVNIIFAVKSTQGFLNGSEVSELLRNLSVVEFSFYLGYPVLQIAEPFQYPQLNLSQLLKSSWVRTVLLGVMEKQLQNEVFQAEMERKLAQLLSEVSTRRRMWRRATVAAGNSVVQVVNVSRLEGDDNPVQLIYFVEDQDGERLSAVKSSDLINKMDLQRAAIILGYRIQGVIAQPVDRVKRPSPESQSNNLWVIVGVVIPVLVVMVIVVILYWKLCRTDKLDFQPDTVANIQQRQKLQIPSVKGFDFAKQHLGQHNKDDILIIHEPAPLPGPLKDHTTPSENGDVPSPKSKIPSKNVRHRGRVSPSDADSTVSEESSERDAGDKTPGAVNDGRSHRAPQSGPPLPSSGNEQHSSASIFEHVDRISRPPEASRRVPSKIQLIAMQPIPAPPVQRPSPADRVAESNKINKEIQTALRHKSEIEHHRNKIRLRAKRRGHYEFPVVDDLSSGDTKERHRVYRRAQMQIDKILDPTASVPSVFIEPRKSSRIKRSPKPRRKHQVNGCPADAEKDRLITTDSDGTYRRPPGVHNSAYIGCPSDPDLPADVQTPSSVELGRYPALPFPASQYIPPQPSIEEARQTMHSLLDDAFALVAPSSQPASTAGVGPGVPPGLPANSTPSQEERRATQWGSFYSPAQTANNPCS",
#                      seq_peptide2 = "DLIRDQGFRGDGGSTTGLSATPPASLPGSLTNVKALQKSPGPQRERKSSSSSEDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATVKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGEFAAFK",
#                      seq_protein = "MPGARRRRRGAAMEGKPRAGVALAPGPSGRRPSARCARRRRPGLLLPGLWLLLLARPASCAPDELSPEQHNLSLYSMELVLKKSTGHSAAQVALTETAPGSQHSSPLHVTAPPSATTFDTAFFNQGKQTKSTADPSIFVATYVSVTSKEVAVNDDEMDNFLPDTHWTTPRMVSPIQYITVSPPGLPREALEPMLTPSLPMVSLQDEEVTSGWQNTTRQPAAYAESASHFHTFRSAFRTSEGIVPTPGRNLVLYPTDAYSHLSSRTLPEIVASLTEGVETTLFLSSRSLMPQPLGDGITIPLPSLGEVSQPPEEVWATSADRYTDVTTVLSQSLEETISPRTYPTVTASHAALAFSRTHSPLLSTPLAFASSASPTDVSSNPFLPSDSSKTSELHSNSALPGPVDNTHILSPVSSFRPYTWCAACTVPSPQQVLATSLMEKDVGSGDGAETLCMTVLEESSISLMSSVVADFSEFEEDPQVFNTLFPSRPIVPLSSRSMEISETSVGISAEVDMSSVTTTQVPPAHGRLSVPASLDPTAGSLSVAETQVTPSSVTTAFFSVITSILLDSSFSVIANKNTPSLAVRDPSVFTPYSLVPSVESSLFSDQERSSFSEHKPRGALDFASSFFSTPPLELSGSISSPSEAPASLSLMPSDLSPFTSQSFSPLVETFTLFDSSDLQSSQLSLPSSTNLEFSQLQPSSELPLNTIMLLPSRSEVSPWSSFPSDSLEFVEASTVSLTDSEAHFTSAFIETTSYLESSLISHESAVTALVPPGSESFDILTAGIQATSPLTTVHTTPILTESSLFSTLTPPDDQISALDGHVSVLASFSKAIPTGTVLITDAYLPSGSSFVSEATPFPLPTELTVVGPSLTPTEVPLNTSTEVSTTSTGAATGGPLDSTLMGDAASQSPPESSAAPPLPSLRPVTAFTLEATVDTPTLATAKPPYVCDITVPDAYLITTVLARRAVQEYIITAIKEVLRIHFNRAVELKVYELFTDFTFLVTSGPFVYTAISVINVLINSKLVRDQTPLILSVKPSFLVPESRFQVQTVLQFVPPSVDTGFCNFTQRIEKGLMTALFEVRKHHQGTYNLTVQILNITISSSRVTPRRGPVNIIFAVKSTQGFLNGSEVSELLRNLSVVEFSFYLGYPVLQIAEPFQYPQLNLSQLLKSSWVRTVLLGVMEKQLQNEVFQAEMERKLAQLLSEVSTRRRMWRRATVAAGNSVVQVVNVSRLEGDDNPVQLIYFVEDQDGERLSAVKSSDLINKMDLQRAAIILGYRIQGVIAQPVDRVKRPSPESQSNNLWVIVGVVIPVLVVMVIVVILYWKLCRTDKLDFQPDTVANIQQRQKLQIPSVKGFDFAKQHLGQHNKDDILIIHEPAPLPGPLKDHTTPSENGDVPSPKSKIPSKNVRHRGRVSPSDADSTVSEESSERDAGDKTPGAVNDGRSHRAPQSGPPLPSSGNEQHSSASIFEHVDRISRPPEASRRVPSKIQLIAMQPIPAPPVQRPSPADRVAESNKINKEIQTALRHKSEIEHHRNKIRLRAKRRGHYEFPVVDDLSSGDTKERHRVYRRAQMQIDKILDPTASVPSVFIEPRKSSRIKRSPKPRRKHQVNGCPADAEKDRLITTDSDGTYRRPPGVHNSAYIGCPSDPDLPADVQTPSSVELGRYPALPFPASQYIPPQPSIEEARQTMHSLLDDAFALVAPSSQPASTAGVGPGVPPGLPANSTPSQEERRATQWGSFYSPAQTANNPCSDLIRDQGFRGDGGSTTGLSATPPASLPGSLTNVKALQKSPGPQRERKSSSSSEDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATVKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGEFAAFK",
#                      path_fasta_tmp = "/media/16TBDisk/Daniela/Fusiones/Fasta-tmp",
#                      path_pfam = "/media/16TBDisk/Daniela/Pfam",
#                      bp1 = 138861139L,
#                      bp2 = 140787584L
#                      )

GetDomains_AllFusions <- function(SeqProt_Transcript) {

  list_domains <- list()
  for(i in 1:nrow(SeqProt_Transcript)) {

    if(SeqProt_Transcript$SeqProtein[i] != "-") {


      domains <- GetDomains(gene1 = SeqProt_Transcript$Gene1[i],
                            gene2 = SeqProt_Transcript$Gene2[i],
                            transcript_id1 = SeqProt_Transcript$transcript_id1[i],
                            transcript_id2 = SeqProt_Transcript$transcript_id2[i],
                            seq_peptide1 = SeqProt_Transcript$Seq_peptide1[i],
                            seq_peptide2 = SeqProt_Transcript$Seq_peptide2[i],
                            #pedazo_peptide1 = SeqProt_Transcript$Pedazo_peptide1[i],
                            #pedazo_peptide2 = SeqProt_Transcript$Pedazo_peptide2[i],
                            seq_protein = SeqProt_Transcript$SeqProtein[i],
                            path_fasta_tmp = "/media/16TBDisk/Daniela/Fusiones/Fasta-tmp",
                            path_pfam = "/media/16TBDisk/Daniela/Pfam",
                            bp1 = SeqProt_Transcript$BP1[i],
                            bp2 = SeqProt_Transcript$BP2[i])

      list_domains[[i]] <- domains
    }
  }
  final_df <- do.call(rbind, list_domains)
  final_df_TK <- final_df[final_df$Domain == "PK_Tyr_Ser-Thr",]
  colnames(final_df_TK)
  final_df_TK <- final_df_TK[, c(1,2,13,14,11,5,6,3,4,7,8,9,10,12)]
  return(list(final_df, final_df_TK))
}



GetDomains <- function(gene1, gene2, transcript_id1, transcript_id2, bp1, bp2, seq_peptide1, seq_peptide2, seq_protein, path_fasta_tmp, path_pfam) {

  Proteins_DB <- data.frame(Transcript_id = character(),Gene = character(), Sequence = character(), stringsAsFactors = FALSE)

  # Añadir la primera fila con gene1 y sequence1
  Proteins_DB <- rbind(Proteins_DB, data.frame(Transcript_id = transcript_id1, Gene = gene1, Sequence = seq_peptide1,  stringsAsFactors = FALSE))
  # Añadir la segunda fila con gene2 y sequence2
  Proteins_DB <- rbind(Proteins_DB, data.frame(Transcript_id = transcript_id2, Gene = gene2, Sequence = seq_peptide2, stringsAsFactors = FALSE))
  Proteins_DB <- unique(Proteins_DB)


  #Escribo el fasta con los pedazos de cada gen que forman las proteínas:

  fasta_file <- file(sprintf("%s/fasta_seq_proteins", path_fasta_tmp), open = "w")
  for (i in 1:nrow(Proteins_DB)) {
    # Escribir el Transcript_id y el Gene en la línea con el símbolo ">"
    writeLines(paste0(">", Proteins_DB$Transcript_id[i], "-", Proteins_DB$Gene[i]), fasta_file)

    # Escribir la secuencia correspondiente en la línea siguiente
    writeLines(Proteins_DB$Sequence[i], fasta_file)
  }

  close(fasta_file)

  # Scan domains with HMMR 3
  pfam_path <- downloadPfam(omicsdo_sof_path = path_pfam)
  system(sprintf("hmmscan --domtblout %s/found-domains.tab %s %s/fasta_seq_proteins", path_fasta_tmp, pfam_path, path_fasta_tmp))
  system(sprintf("cat %s/found-domains.tab | grep -v '^#' | sed 's/  */\t/g' | cut -f 1,2,4,20,21 > %s/found-domains-extract.tab", path_fasta_tmp, path_fasta_tmp))

  domains_df <- read.table(sprintf("%s/found-domains-extract.tab", path_fasta_tmp), header = FALSE, sep = "\t",
                           col.names = c("Domain", "DomainID", "Gen", "Start", "End"))

  if(nrow(domains_df) == 0) {
    domains_df[1, ] <- "-"
    message(sprintf("NO domains found for fusion:
                    gene1 = %s,
                    gene2 = %s,
                    transcript_id1 = %s,
                    transcript_id2 = %s,
                    bp1 = %s,
                    bp2 = %s", gene1, gene2, transcript_id1, transcript_id2, bp1, bp2))
    domains_df_1 <- data.frame("Gen" = gene1,
                               "TranscriptID" = transcript_id1,
                               "BP" = bp1,
                               "Seq_protein" = seq_protein,
                               "Domain" = "-",
                               "DomainID" = "-",
                               "Start" = "-",
                               "End" = "-",
                               "seq_peptide" = seq_peptide1,
                               "Domain_Sequence" = "-",
                               "Perc_Conserved_Domain" =  0,
                               "Conserved_Domain" = "-",
                               "Coords_Conserved_Domain" = "-",
                               "Residues_CBDOCK2" = "-")
    domains_df_2 <- data.frame("Gen" = gene2,
                               "TranscriptID" = transcript_id2,
                               "BP" = bp2,
                               "Seq_protein" = seq_protein,
                               "Domain" = "-",
                               "DomainID" = "-",
                               "Start" = "-",
                               "End" = "-",
                               "seq_peptide" = seq_peptide2,
                               "Domain_Sequence" = "-",
                               "Perc_Conserved_Domain" =  0,
                               "Conserved_Domain" = "-",
                               "Coords_Conserved_Domain" = "-",
                               "Residues_CBDOCK2" = "-")

    domains_df <- rbind(domains_df_1, domains_df_2)

  } else {

    domains_df$TranscriptID <- sapply(strsplit(domains_df$Gen, split = "-"), `[`, 1)
    domains_df$Gen <- sapply(strsplit(domains_df$Gen, split = "-"), `[`, 2)
    colnames(domains_df)
    domains_df <- domains_df[, c(3,6,1,2,4,5)]

    domains_df$seq_peptide <- "-"
    domains_df$Domain_Sequence <- "-"
    domains_df$BP <- "-"
    domains_df$Seq_protein <- seq_protein
    i=1
    for(i in 1:nrow(domains_df)) {
      if(domains_df$TranscriptID[i] == transcript_id1) {
        domains_df$seq_peptide[i] <- seq_peptide1
        domains_df$BP[i] <- bp1
      } else if(domains_df$TranscriptID[i] %in% transcript_id2) {
        domains_df$seq_peptide[i] <- seq_peptide2
        domains_df$BP[i] <- bp2
      }
      domains_df$Domain_Sequence[i] <- substr(domains_df$seq_peptide[i], domains_df$Start[i], domains_df$End[i])
    }

    #Agregar coordenadas:
    #write.xlsx(domains_df, file = sprintf("%s/Domains_df.xlsx", path_dir))

    #Poner el dominio en el fusions_report: --------------------------------------
    i=1
    domains_df$Perc_Conserved_Domain <- "-"
    domains_df$Conserved_Domain <- "-"
    domains_df$Coords_Conserved_Domain <- "-"
    domains_df$Residues_CBDOCK2 <- "-"

    for( i in 1:nrow(domains_df)) {

      #ANALIZO EL DOMINIO DEL GEN 1:
      domain <- domains_df$Domain_Sequence[i]

      message("Analizando el porcentaje de conservación del dominio")
      library(Biostrings)
      conserved_domain <- get_longest_common_substring(domain, seq_protein)
      porcentaje_conservacion <- round(nchar(conserved_domain) / nchar(domain) * 100)

      domains_df$Perc_Conserved_Domain[i] <- porcentaje_conservacion
      domains_df$Conserved_Domain[i] <- conserved_domain

      # Encontrar coordenadas del conserved_dominio dentro de la proteina:
      start_domain <- regexpr(conserved_domain, seq_protein)[1]
      end_domain <- start_domain + nchar(conserved_domain) -1
      domains_df$Coords_Conserved_Domain[i] <- paste0(start_domain, "-", end_domain)

      # Generar residuos para CBDOCK2:
      positions <- seq(start_domain, end_domain)
      amino_acids <- strsplit(domain, "")[[1]]
      domains_df$Residues_CBDOCK2[i] <- paste0(positions, ":", amino_acids, collapse = ",")
    }

    colnames(domains_df)
    domains_df <- domains_df[, c(1,2,9,10,3,4,5,6,7,8,11,12,13,14)]

  }
  #write.xlsx(domains_df, file = sprintf("%s/Fusions_Annot_Domains.xlsx", path_dir))
  return(domains_df)
}

#ANOTACION DE FUSIONES -----------------------------------------------------------------------------------

#' @title AnotateFusions
#' @description Generates
#' @param fusions_report df
#' @param protein_df df
#' @export
#' @import readxl
#' @import readr

#gene1 = "KIAA1549"
#gene2 = "BRAF"
#bp1 = 138867975L
#bp2 = 140787584L
#transcript_id1 ="ENST00000422774"
#transcript_id2 = "ENST00000496384"
#library(readr)
#BED_Anotador <- read_rds("/media/16TBDisk/Daniela/OMICsdo/data/BED_Anotador.RDS")
#protein_df <- read_rds("/media/16TBDisk/Daniela/OMICsdo/data/DB_Anotador_Prot_DNA.RDS")


#Encuentra la posición del bp dentro de la secuencia del DNA: ----
cdna_pos_from_bp <- function(exons_df, transcript_id, bp, dna_seq) {
  ex <- subset(exons_df, transcript == transcript_id)
  ex$exonNumber <- as.numeric(ex$exonNumber)

  #if(unique(ex$strand == "+")) {
  ex <- ex[order(ex$exonNumber), ]
  #} else {
  #  ex <- ex[order(ex$exonNumber, decreasing = TRUE), ]
  #}

  lens <- ex$end - ex$start + 1L

  # Encontrar el exón que contiene el bp
  hit <- which(bp >= ex$start & bp <= ex$end)
  if (length(hit) != 1) stop("El bp no cae dentro de un único exón de ese transcripto.")

  # Acumulado previo
  cum_prev <- if (hit == 1) 0L else sum(lens[1:(hit-1)])

  if (unique(ex$strand) == "+") {
    offset <- (bp - ex$start[hit]) + 1L
  } else {
    offset <- (ex$end[hit] - bp) + 1L
  }

  pos <- cum_prev + offset
  if (pos < 1L || pos > nchar(dna_seq)) warning("La posición calculada excede la longitud de dna_seq.")
  return(pos)
}
#--------


#Traduce de DNA a Peptidos respetando codon de start y de stop:
from_dna_to_peptide <- function(lado_5prime) {

  # 1. Buscar primer ATG
  start_pos <- as.integer(start(matchPattern("ATG", lado_5prime))[1])

  if (!is.na(start_pos)) {
    # Recortar desde el primer ATG
    seq_from_start <- subseq(lado_5prime, start = start_pos)

    # 2. Buscar stop codones en el mismo marco
    stop_codons <- c("TAA", "TAG", "TGA")
    stop_matches <- lapply(stop_codons, function(stop) {
      start(matchPattern(stop, seq_from_start))
    })
    stop_positions <- sort(unlist(stop_matches))

    # Filtrar solo los que están en el mismo marco que ATG
    stop_in_frame <- stop_positions[(stop_positions %% 3) == 1]

    if (length(stop_in_frame) > 0) {
      # Tomar el primer stop en marco
      stop_pos <- stop_in_frame[1] + 2  # incluir las 3 bases del stop

      # 3. Recortar desde ATG hasta STOP
      coding_seq <- subseq(seq_from_start, start = 1, end = stop_pos)

      # 4. Traducir a aminoácidos
      peptide <- translate(DNAString(coding_seq))
      peptide
    } else {
      message("No se encontró codón STOP en el mismo marco.")
      peptide <- translate(DNAString(seq_from_start))
    }
  } else {
    message("No se encontró codón de inicio ATG.")
    peptide <- NA
  }
  return(peptide)
}
#-------------------------------
from_dna_to_peptide2 <- function(lado_5prime) {

  # 1. Buscar primer ATG
  #start_pos <- as.integer(start(matchPattern("ATG", lado_5prime))[1])

  #if (!is.na(start_pos)) {
  # Recortar desde el primer ATG
  #seq_from_start <- subseq(lado_5prime, start = start_pos)
  seq_from_start <- lado_5prime

  # 2. Buscar stop codones en el mismo marco
  stop_codons <- c("TAA", "TAG", "TGA")
  stop_matches <- lapply(stop_codons, function(stop) {
    start(matchPattern(stop, seq_from_start))
  })
  stop_positions <- sort(unlist(stop_matches))

  # Filtrar solo los que están en el mismo marco que ATG
  stop_in_frame <- stop_positions[(stop_positions %% 3) == 1]

  if (length(stop_in_frame) > 0) {
    # Tomar el primer stop en marco
    stop_pos <- stop_in_frame[1] + 2  # incluir las 3 bases del stop

    # 3. Recortar desde ATG hasta STOP
    coding_seq <- subseq(seq_from_start, start = 1, end = stop_pos)

    # 4. Traducir a aminoácidos
    peptide <- translate(DNAString(coding_seq))
    peptide
  } else {
    message("No se encontró codón STOP en el mismo marco.")
    peptide <- translate(DNAString(seq_from_start))
  }
  #} else {
  #  message("No se encontró codón de inicio ATG.")
  #  peptide <- NA
  #}
  return(peptide)
}
#-------------------------------
#bp1 = 138861139L
#bp2 = 140787584L
library(Biostrings)



Anotate_All_Transcripts <- function(gene1, gene2, bp1, bp2, protein_df) {
  transcripts_gene1 <- protein_df$Transcript_clean[which(protein_df$gene_name == gene1)]
  transcripts_gene1 <- transcripts_gene1[!duplicated(transcripts_gene1)]

  transcripts_gene2 <- protein_df$Transcript_clean[which(protein_df$gene_name == gene2)]
  transcripts_gene2 <- transcripts_gene2[!duplicated(transcripts_gene2)]

  # Generar todas las combinaciones posibles
  combinaciones <- expand.grid(transcript_id1 = transcripts_gene1,
                               transcript_id2 = transcripts_gene2,
                               stringsAsFactors = FALSE)

  #comb = 1
  SeqProt_Transcript <- combinaciones
  SeqProt_Transcript$Gene1 <- gene1
  SeqProt_Transcript$Gene2 <- gene2
  SeqProt_Transcript$BP1 <- bp1
  SeqProt_Transcript$BP2 <- bp2

  SeqProt_Transcript$SeqProtein <- "-"
  SeqProt_Transcript$Pedazo_peptide1 <- "-"
  SeqProt_Transcript$Pedazo_peptide2 <- "-"

  SeqProt_Transcript$Seq_peptide1 <- "-"
  SeqProt_Transcript$Seq_peptide2 <- "-"
  SeqProt_Transcript$n_AA <- "-"
  SeqProt_Transcript$Seq_dna1 <- "-"
  SeqProt_Transcript$Seq_dna2 <- "-"
  SeqProt_Transcript$Comment <- "-"

  BED_Anotador <- read_rds("/media/16TBDisk/Daniela/OMICsdo/data/BED_Anotador.RDS")
  #protein_df   <- read_rds("/media/16TBDisk/Daniela/OMICsdo/data/DB_Anotador_Prot_DNA.RDS")

  comb=1
  for(comb in 1:nrow(combinaciones)) {
    t1 <- combinaciones$transcript_id1[comb]
    t2 <- combinaciones$transcript_id2[comb]

    tryCatch(
      {

        o <- AnotateFusions(
          gene1          = gene1,
          gene2          = gene2,
          bp1            = bp1,
          bp2            = bp2,
          transcript_id1 = t1,
          transcript_id2 = t2,
          BED_Anotador   = BED_Anotador,
          protein_df     = protein_df
        )

        SeqProt_Transcript$SeqProtein[comb] <- o[[1]]
        SeqProt_Transcript$Pedazo_peptide1[comb] <- o[[2]] # pedazo_seq 1 completo
        SeqProt_Transcript$Pedazo_peptide2[comb] <-  o[[3]] # pedazo seq 2 completo
        SeqProt_Transcript$Seq_peptide1[comb] <- o[[4]]
        SeqProt_Transcript$Seq_peptide2[comb] <-  o[[5]]
        SeqProt_Transcript$n_AA[comb] <- nchar(o[[1]])
        SeqProt_Transcript$Seq_dna1[comb] <-  o[[6]]
        SeqProt_Transcript$Seq_dna2[comb] <-  o[[7]]

      },
      error = function(e) {
        message("Error en AnotateFusions: ", e$message)
        error <- e$message
      }
    )

  }

  filas_errores <- which(SeqProt_Transcript$SeqProtein == "-")
  SeqProt_Transcript$Comment[filas_errores] <- "El bp no cae dentro de un exón de ese transcripto."
  colnames(SeqProt_Transcript)
  SeqProt_Transcript <- SeqProt_Transcript[, c(3,1,5,4,2,6,7,12, 8, 9,10,11,13, 14, 15)]
}


#Ejemplo del FusionGDB
#o <- AnotateFusions(gene1 = "KIAA1549",
#                    gene2 = "BRAF",
#                    bp1 = 138545884L,
#                    bp2 = 140487384L,
#                    transcript_id1 ="ENST00000440172",
#                    transcript_id2 = "ENST00000288602",
#                    BED_Anotador <- read_rds("/media/16TBDisk/Daniela/OMICsdo/data/BED_Anotador.RDS"),
#                    protein_df <- read_rds("/media/16TBDisk/Daniela/OMICsdo/data/DB_Anotador_Prot_DNA.RDS")
#)

#seq_total <- o[[1]]
#seq1 <- o[[2]]
#seq2 <- o[[3]]

#seq_total == seq_prot_fusionada

AnotateFusions <- function(gene1, gene2, transcript_id1, transcript_id2, bp1, bp2, BED_Anotador, protein_df) {

  library(stringr)

  #HAGO ANALISIS POR TRANSCRIPT ID: ---------------------------------------------
  if( !is.na(transcript_id1)) {

    if(transcript_id1 %in% protein_df$Transcript | transcript_id1 %in% protein_df$Transcript_clean) {

      transcript_sequence1 <- protein_df$Sequence[which(protein_df$Transcript == transcript_id1 | protein_df$Transcript_clean == transcript_id1 )]
      if(any(duplicated(transcript_sequence1))){
        transcript_sequence1 <- transcript_sequence1[!duplicated(transcript_sequence1)]
      }

      dna_sequence1 <- protein_df$DNA_Sequence[which(protein_df$Transcript == transcript_id1 | protein_df$Transcript_clean == transcript_id1 )]
      if(any(duplicated(dna_sequence1))){
        dna_sequence1 <- dna_sequence1[!duplicated(dna_sequence1)]
      }

      #bed1 <- protein_df[protein_df$gene_name == gene1 & (protein_df$Transcript == transcript_id1 | protein_df$Transcript_clean == transcript_id1),]
      #bed1 <- bed1[-which(is.na(bed1$Transcript_clean)),]

      exons1 <- BED_Anotador[BED_Anotador$geneName == gene1 & BED_Anotador$type == "exon" & BED_Anotador$transcript == transcript_id1,]

      #Encuentro la posición del bp en la secuencia del DNA:
      pos_cdna <- cdna_pos_from_bp(exons_df = exons1,
                                   transcript_id = transcript_id1,
                                   bp = bp1,
                                   dna_seq = dna_sequence1)

      #Armo pedazo de seq del DNA del gen 1:
      lado_5prime <- subseq(dna_sequence1, 1, pos_cdna) # hasta el bp (incluye bp)
      #lado_5prime <- subseq(dna_sequence1, 1, pos_cdna-1) # hasta el bp (incluye bp)

      nchar(lado_5prime)
      nchar(dna_sequence1)

      #Armo pedazo de seq de peptidos del gen 1:
      lado_5prime_aa <- from_dna_to_peptide(lado_5prime = DNAString(lado_5prime))
      pedazo_seq1_completo <- as.character(lado_5prime_aa)
      nchar(pedazo_seq1_completo)
      nchar(transcript_sequence1)
      pedazo_seq1_completo <- gsub("\\*", "", pedazo_seq1_completo)


    } else {
      message("NO se encuentra el transcript id1 ni en el clean")
      pedazo_seq1_completo <- "."
    }
  }

  #lo mismo para el trasncript id 2 ----------------------------------------
  if( !is.na(transcript_id2)) {

    if(transcript_id2 %in% protein_df$Transcript_clean) {

      transcript_sequence2 <- protein_df$Sequence[which(protein_df$Transcript_clean == transcript_id2 )]
      if(any(duplicated(transcript_sequence2))){
        transcript_sequence2 <- transcript_sequence2[!duplicated(transcript_sequence2)]
      }

      dna_sequence2 <- protein_df$DNA_Sequence[which(protein_df$Transcript_clean == transcript_id2 )]
      if(any(duplicated(dna_sequence2))){
        dna_sequence2 <- dna_sequence2[!duplicated(dna_sequence2)]
      }

      exons2 <- BED_Anotador[BED_Anotador$geneName == gene2 & BED_Anotador$type == "exon"& BED_Anotador$transcript == transcript_id2,]
      #exons2$dif <- exons2$end - exons2$start
      #length_exons <- sum(exons2$dif)
      #nchar(dna_sequence2)

      #Prueba:----
      #dna_aa <- from_dna_to_peptide2(lado_5prime = DNAString(dna_sequence2))
      #dna_aa <- as.character(dna_aa)
      #nchar(dna_aa)

      #-------------

      #Encuentro la posición del bp en la secuencia del DNA:
      pos_cdna2 <- cdna_pos_from_bp(exons_df = exons2,
                                    transcript_id = transcript_id2,
                                    bp = bp2,
                                    dna_seq = dna_sequence2)

      #Armo pedazo de seq del DNA del gen 2:
      #lado_3prime <- subseq(dna_sequence2, pos_cdna + 1, nchar(dna_sequence2)) # después del bp
      lado_3prime <- subseq(dna_sequence2, pos_cdna2, nchar(dna_sequence2)) # después del bp

      #con esto anda bien:
      #lado_3prime <- subseq(dna_sequence2, pos_cdna -107, nchar(dna_sequence2)) # después del bp

      #Armo pedazo de seq de peptidos del gen 2:
      lado_3prime_aa <- from_dna_to_peptide2(lado_5prime = DNAString(lado_3prime))
      pedazo_seq2_completo <- as.character(lado_3prime_aa)
      #pedazo_seq2_completo
      pedazo_seq2_completo <- gsub("\\*", "", pedazo_seq2_completo)

    } else {
      message("NO se encuentra el transcript id2 ni en el clean")
      pedazo_seq2_completo <- "."
    }
  }

  seq_prot_fusionada <- paste0(pedazo_seq1_completo, pedazo_seq2_completo)
  #seq_prot_fusionada <- toupper(seq_prot_fusionada)

  return(list(seq_prot_fusionada, pedazo_seq1_completo, pedazo_seq2_completo, transcript_sequence1, transcript_sequence2, dna_sequence1, dna_sequence2))
}


