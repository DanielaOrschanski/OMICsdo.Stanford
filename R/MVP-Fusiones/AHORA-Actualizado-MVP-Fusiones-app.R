# FLENI Fusiones - Shiny App para Visualización y Anotación de Fusiones

























# Cargar librerías necesarias
#install.packages("echarts4r")

library(shiny)
library(shinydashboard)
library(DT)
library(readxl)
library(dplyr)
library(shinyFiles)
library(shinycssloaders)
library(plotly)
library(echarts4r)

# Usar renderEchartsforR para que los graficos sean mas esteticos

# Ruta al archivo
ruta_archivo <- "/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/Act_ReportFusions_TodasCohortes.xlsx"
datos <- read_excel(ruta_archivo)
comas <- datos[which(grepl(",", datos$gene1) | grepl(",", datos$gene2)),]
# Ruta del archivo Excel
ruta_excel <- "/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/Act_ReportFusions_TodasCohortes.xlsx"
ruta_excel_genes <- "/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/ReportGenes_TodasCohortes.xlsx"

datos_genes <-read.xlsx(ruta_excel_genes)

library(readr)
protein_df <- read_rds("/media/16TBDisk/Daniela/OMICsdo/data/DB_Anotador_Prot_DNA.RDS")
load("/media/16TBDisk/Daniela/OMICsdo/Grafico-Fusiones-Propio/exons-arriba.RData")

Metadata_FLENI <- read_excel("/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/TablaPacientesRNAseq_Fusiones.xlsx")
Metadata_FLENI$Cohort <- "FLENI"
colnames(Metadata_FLENI)[1] <- "ID"

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Fusion.AR-T"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("by Patient", tabName = "fusiones", icon = icon("dna")),
      menuItem("by Fusion", tabName = "analisis", icon = icon("microscope")),
      menuItem("by Gene", tabName = "genes", icon = icon("glyphicon glyphicon-asterisk")),
      menuItem("QC", tabName = "qc", icon = icon("chart-bar"))
    )
  ),
  dashboardBody(
    tabItems(
      # --------------------- Dashboard ----------------------------------------------
      tabItem(tabName = "dashboard",

              fluidRow(
                valueBoxOutput("num_muestras", width = 3),
                valueBoxOutput("num_fusiones", width = 3),
                valueBoxOutput("num_cohortes", width = 3),
                valueBoxOutput("num_tipos_cancer", width = 3)
              ),

              fluidRow(
                box(width = 12, title = "Samples' Distributions", status = "primary", solidHeader = TRUE,
                    fluidRow(
                      column(width = 6,
                             plotlyOutput("graficoCohorte")
                      ),
                      column(width = 6,
                             plotlyOutput("graficoCancer")
                      )
                    )
                  )
                ),

              fluidRow(
                box(width = 12, title = "Fusions' Distributions", status = "primary", solidHeader = TRUE,
                    # Primera fila de gráficos
                    fluidRow(
                      column(width = 6,
                             plotlyOutput("FusCohorte")
                      ),
                      column(width = 6,
                             plotlyOutput("FusCancer")
                      )
                    ),

                )
              ),

              # Gŕaficos dinamicos ---------------------------------------
              fluidRow(
                box(width = 12, title = "Fusions' Explorer", status = "primary", solidHeader = TRUE,
                  fluidRow(
                    column(
                      width = 4,
                      selectInput(
                        inputId = "columnasAgrupar",
                        label = "Selecciona columnas para agrupar:",
                        choices = names(read_excel(ruta_excel)), # todas las columnas del excel
                        multiple = TRUE,
                        selected = "Cohort" # valor inicial
                      )
                    ),
                    column(4,
                           selectInput("columnaFiltrarGen", "Selecciona columna de gen:",
                                       choices = c("gene1", "gene2"))
                    ),
                    column(4,
                           uiOutput("seleccionGen")  # selectInput dinámico de genes
                    ),

                    #Definir eje Y: fusiones por ID o fusiones totales --------
                    column(
                      width = 4,
                      checkboxInput(
                        inputId = "usarTotalFusiones",
                        label = "Nº total de fusiones (barra) en vez de Fusiones por ID (boxplot)",
                        value = FALSE
                      )
                    ),

                    column(width = 12, plotlyOutput("BoxplotCancer")),
                  )
                )
              ),


              fluidRow(
                box(width = 12, title = "TK' Distributions", status = "primary", solidHeader = TRUE,
                    fluidRow(
                      column(width = 6,
                             plotlyOutput("TK_Cohorte")
                      ),
                      column(width = 6,
                             plotlyOutput("TK_100_Cohorte")
                      )
                    )
                )
              ),

              #Tabla con todas las fusiones
              fluidRow(
                box(width = 12,
                    title = "All Fusions",
                    status = "primary",
                    solidHeader = TRUE,
                    DTOutput("tabla")
                )
              ),
            ),
    # ---------------- by Patient ---------------- -----------------

    tabItem(tabName = "fusiones",
            fluidRow(
              box(width = 12, title = "Patient Explorer", status = "primary", solidHeader = TRUE,

                  # Selectores dinámicos
                  fluidRow(
                    column(width = 4,
                           selectInput("selectCohorte", "Select Cohort:", choices = NULL)
                    ),
                    column(width = 4,
                           selectInput("selectID", "Select Patient ID:", choices = NULL)
                    )
                  ),
                  fluidRow(
                    column(width = 4,
                           checkboxInput("fusionesConocidas", "Apply gene panel (known fusions)", value = FALSE)
                    )
                  ),
                  fluidRow(
                    column(width = 4,
                           checkboxInput("fusionesTK", "Keep Tirosin Kinase", value = FALSE)
                    )
                  ),

                  # Tabla de fusiones filtradas
                  DTOutput("tablaPaciente"),
              )
            ),

            # Acá va lo dinámico: una UI para cada fusión seleccionada
            fluidRow(
              box(width = 12, title = "Fusion details", status = "info", solidHeader = TRUE,
                  uiOutput("fusionPanels")
              )
            ),

            # Metadata by patient --------------
            fluidRow(
              box(width = 12, title = "Clinical Information", status = "info", solidHeader = TRUE,
                  DTOutput("tablaPacienteMeta")
              )
            ),

            # QC by patient --------------
            fluidRow(
              box(width = 6, title = "QC Information", status = "info", solidHeader = TRUE,
                  DTOutput("tablaPacienteQC")
              ),
              box(width = 6, title = "QC Percentiles", status = "info", solidHeader = TRUE,
                   plotlyOutput("QCBarras"),
                  br(),
                  htmlOutput("explicacion_score_qc")
              )
            ),

            # 🔹 Botón de descarga al final de toda la pestaña
            fluidRow(
              box(width = 12, title = NULL, solidHeader = FALSE, status = "primary",
                  downloadButton("downloadReport", "Download Report", class = "btn-primary")
              )
            )

        ),
    #----------------------------------------------
    # by Fusion
    tabItem(tabName = "analisis",
            fluidRow(
              column(6,
                     selectInput("gen5", "Gen 5'",
                                 choices = unique(protein_df$gene_name),
                                 selected = NULL)
              ),
              column(6,
                     selectInput("gen3", "Gen 3'",
                                 choices = unique(protein_df$gene_name),
                                 selected = NULL)
              )
            ),
            fluidRow(
              column(6,
                     selectizeInput("bp1", "Breakpoint 1",
                                    choices = NULL,  # Se completa dinámicamente
                                    options = list(create = TRUE))  # permite ingresar otro número
              ),
              column(6,
                     selectizeInput("bp2", "Breakpoint 2",
                                    choices = NULL,
                                    options = list(create = TRUE))
              )
            ),

            br(),
            # 👇 Aquí agregamos el botón
            fluidRow(
              column(12,
                     actionButton("predict_btn", "Predict Protein",
                                  icon = icon("dna"), class = "btn-primary")
              )
            ),

            br(),
            fluidRow(
              column(12, h3("Annotated Transcripts")),
              column(12, DT::dataTableOutput("tabla_fusion_transcript"))
            ),
            fluidRow(
              column(12, h3("Domains")),
              column(12, DT::dataTableOutput("tabla_fusion_domains"))
            )
        ),

    # QC -----------------------------------------
    tabItem(tabName = "qc",

            # Selectores dinámicos
            fluidRow(
              column(width = 4,
                     selectInput("selectCohorteQC", "Select Cohort:", choices = NULL)
              ),
              column(width = 4,
                     selectInput("selectIDQC", "Select Patient ID:", choices = NULL))
            ),

            # --- Graficos ---
            fluidRow(
              box(title = "Métricas QC Secuenciación", width = 12, status = "primary", solidHeader = TRUE,
                  plotlyOutput("qc_secuenciacion_plot", height = "400px"))
            ),
            fluidRow(
              box(title = "Métricas QC Alineamiento", width = 12, status = "success", solidHeader = TRUE,
                  plotlyOutput("qc_alineamiento_plot", height = "400px"))
            )
    ),
    # ---------------- by Gene ---------------- -----------------

    tabItem(tabName = "genes",
            fluidRow(
              box(width = 12, title = "Top Fused Genes", status = "primary", solidHeader = TRUE,

                  # Selectores dinámicos
                  fluidRow(
                    column(width = 4,
                           selectInput("selectCohorteGenes", "Select Cohort:", choices = NULL)
                    )
                  ),
                  fluidRow(
                    column(width = 4,
                           checkboxInput("genesConocidos", "Apply gene panel (known fusion genes)", value = FALSE)
                    )
                  ),

                  # Tabla de fusiones filtradas
                  DTOutput("tablaGenes"),
              )
            )

      )
    #----------------------------------------------

#------------------------------------------
    )
  )
)

#FUNCIONES A USAR EN EL SERVER --------------------------

# PARA BY PATIENT ################################

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

generate_PedazoSeqCompleto <- function(pedazo_seq, seq, gen) {
  if(is.na(pedazo_seq)) {
    message("pedazo es NA")
    pedazo_seq <- ""
  }
  if(!grepl(pedazo_seq, seq)) {
    message("Convierto minusuclas")
    pedazo_seq <- toupper(pedazo_seq)
  }
  if(grepl("\\?", pedazo_seq)) {
    pedazo_seq <- gsub("\\?", "", pedazo_seq)
  }

  if(grepl("\\*", pedazo_seq)) {
    message("Se encontró un codón de STOP")
    pedazo_seq <- gsub("\\*", "", pedazo_seq)
    pedazo_seq_completo <- ""
  } else {
    if(!grepl(pedazo_seq, seq)) {
      message("Elimino el ultimo AA")
      pedazo_seq <- substr(pedazo_seq, 1, nchar(pedazo_seq) - 1)
      agregar_AAfinal <- "SI"
    }
    if(!grepl(pedazo_seq, seq)) {
      message("No se encuentra el pedazo de peptide en la secuencia completa")
      pedazo_seq_completo <-  ""
    } else {
      message("Si se encontró el pedazo en la secuencia")
      #Si es gen 1 se queda con la parte de adelante de la secuencia.
      # Si es gen 2 se queda con la parte de atras de la secuencia
      pedazo_seq_completo <- ifelse(gen == 1,
                                    str_split(seq, pedazo_seq)[[1]][1],
                                    str_split(seq, pedazo_seq)[[1]][2]
      )
    }
  }

  if(is.na(pedazo_seq_completo) | pedazo_seq== "") {
    pedazo_seq_completo <- ""
    return("NO SE ENCUENTRA EL PEDAZO DE PEPTIDE EN LA SECUENCIA DEL GEN")
  }

  #Si es gen 1 se pega secuencia + pedazo.
  # Si es gen 2 se pega pedazo + secuencia.
  pedazo_seq_completo <- ifelse(gen == 1,
                                paste0(pedazo_seq_completo, pedazo_seq),
                                paste0(pedazo_seq, pedazo_seq_completo))

  return(pedazo_seq_completo)
}


getDrugs <- function(protein_df, gene) {

  if(gene %in% protein_df$gene_name) {
    drugs_gene <- protein_df[which(protein_df$gene_name == gene),]
    if(all(is.na(drugs_gene$drug_name))) {

      drugs_gene <- drugs_gene[1, which(colnames(drugs_gene) %in% c("gene_name", "drug_name"))]
      drugs_gene$drug_name <- "No drugs approved for this gene"

    } else {
      drugs_gene <- drugs_gene[, which(colnames(drugs_gene) %in% c("drug_name", "interaction_score", "interaction_type", "approved", "immunotherapy", "anti_neoplastic"))]
      if(any(duplicated(drugs_gene))) {
        drugs_gene <- drugs_gene[!duplicated(drugs_gene),]
      }
      drugs_gene <- drugs_gene[which(drugs_gene$approved == "TRUE"),]
      drugs_gene <- drugs_gene[which(drugs_gene$immunotherapy == "FALSE"),]
      drugs_gene <- drugs_gene[, c(3,1,2,4,6)]
      rownames(drugs_gene) <- NULL
      drugs_gene <- drugs_gene[order(drugs_gene$interaction_type),]
    }
  } else {
    drugs_gene <- data.frame(gene_name = gene, drug_name = "No drugs approved for this gene")
  }

  return(drugs_gene)
}


AnotateFusion <- function(protein_df, transcript_id1, transcript_id2, peptide_sequence) {

  library(stringr)
  agregar_AAfinal <- "NO"

  if(transcript_id1 %in% protein_df$Transcript) {
    sequence1 <- protein_df$Sequence[which(protein_df$Transcript == transcript_id1)]
    if(any(duplicated(sequence1))) {
      sequence1 <- sequence1[!duplicated(sequence1)]
    }
    pedazo_seq1 <- str_split(peptide_sequence, "\\|")[[1]][1]
    pedazo_seq1_completo <- generate_PedazoSeqCompleto(pedazo_seq = pedazo_seq1,
                                                       seq = sequence1,
                                                       gen = 1)
  } else { # SI NO ENCONTRAMOS EL EXACTO TRANSCRIPT ID: --------------------------

    sequence1 <- "No se encuentra el transcript id1"
    message("No se encuentra el transcript id1")
    transcript_id1_clean <- strsplit(transcript_id1, split="\\.")[[1]][1]

    if(transcript_id1_clean %in% protein_df$Transcript_clean) {
      message("SI se encuentra el transcript clean")

      nuevo_transcript_id1 <- protein_df$Transcript[which(protein_df$Transcript_clean == transcript_id1_clean)]
      sequence1 <- protein_df$Sequence[which(protein_df$Transcript == nuevo_transcript_id1)]
      if(any(duplicated(sequence1))) {
        sequence1 <- sequence1[!duplicated(sequence1)]
      }
      pedazo_seq1 <- str_split(peptide_sequence, "\\|")[[1]][1]
      pedazo_seq1_completo <- generate_PedazoSeqCompleto(pedazo_seq = pedazo_seq1,
                                                         seq = sequence1,
                                                         gen = 1)
    }
    else {
      message("NO se encuentra el transcript id1 ni en el clean")
      pedazo_seq1_completo <- "."
    }
  }

  #lo mismo para el trasncript id 2 ----------------------------------------

  if(transcript_id2 %in% protein_df$Transcript) {
    sequence2 <- protein_df$Sequence[which(protein_df$Transcript == transcript_id2)]
    if(any(duplicated(sequence2))) {
      sequence2 <- sequence2[!duplicated(sequence2)]
    }
    agregar_AAfinal <- "NO"
    pedazo_seq2 <- str_split(peptide_sequence, "\\|")[[1]][2]
    pedazo_seq2_completo <- generate_PedazoSeqCompleto(pedazo_seq = pedazo_seq2,
                                                       seq = sequence2,
                                                       gen = 2)

  } else {
    sequence2 <- "No se encuentra el transcript id2"
    message("No se encuentra el transcript id2")
    transcript_id2_clean <- strsplit(transcript_id2, split="\\.")[[1]][1]

    if(transcript_id2_clean %in% protein_df$Transcript_clean) {
      message("SI se encuentra el transcript clean")
      nuevo_transcript_id2 <- protein_df$Transcript[which(protein_df$Transcript_clean == transcript_id2_clean)]
      sequence2 <- protein_df$Sequence[which(protein_df$Transcript == nuevo_transcript_id2)]
      if(any(duplicated(sequence2))) {
        sequence2 <- sequence2[!duplicated(sequence2)]
      }
      pedazo_seq2 <- str_split(peptide_sequence, "\\|")[[1]][2]
      pedazo_seq2_completo <- generate_PedazoSeqCompleto(pedazo_seq = pedazo_seq2,
                                                         seq = sequence2,
                                                         gen = 2)
    } else {
      message("NO se encuentra el transcript id2 ni en el clean")
      pedazo_seq2_completo <- "."
    }
  }

  seq_prot_fusionada <- paste0(pedazo_seq1_completo, pedazo_seq2_completo)
  seq_prot_fusionada <- toupper(seq_prot_fusionada)

  return(list(seq_prot_fusionada, sequence1, sequence2, pedazo_seq1_completo, pedazo_seq2_completo))
}

#transcript_id = "ENST00000646891"
#gene = "BRAF"
#sequence_gene = "MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEHIEALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTVTSSSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRDSLKKALMMRGLIPECCAVYRIQDGEKKPIGWDTDISWLTGEELHVEVLENVPLTTHNFVRKTFFTLAFCDFCRKLLFQGFRCQTCGYKFHQRCSTEVPLMCVNYDQLDLLFVSKFFEHHPIPQEEASLAETALTSGSSPSAPASDSIGPQILTSPSPSKSIPIPQPFRPADEDHRNQFGQRDRSSSAPNVHINTIEPVNIDDLIRDQGFRGDGGSTTGLSATPPASLPGSLTNVKALQKSPGPQRERKSSSSSEDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATEKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGAFPVH"
#seq_prot = "MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEHIEALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTVTSSSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRDSLKKALMMRGLIPECCAVYRIQDGEKKPIGWDTDISWLTGEELHVEVLENVPLTTHNFVRKTFFTLAFCDFCRKLLFQGFRCQTCGYKFHQRCSTEVPLMCVNYDQLDLLFVSKFFEHHPIPQEEASLAETALTSGSSPSAPASDSIGPQILTSPSPSKSIPIPQPFRPADEDHRNQFGQRDRSSSAPNVHINTIEPVNIDDLIRDQGFRGDGGSTTGLSATPPASLPGSLTNVKALQKSPGPQRERKSSSSSEDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATEKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGAFPVH"

GetDomain <- function(transcript_id, gene, sequence_gene, seq_pedazo, seq_prot, path_dir = "/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/fasta_tmp") {

  Proteins_DB <- data.frame(Transcript_id = transcript_id, Gene = gene, Sequence = sequence_gene, stringsAsFactors = FALSE)
  if(any(Proteins_DB$Transcript_id == ".")) { Proteins_DB <- Proteins_DB[-which(Proteins_DB$Transcript_id == "."),] }
  if(any(Proteins_DB$Sequence == "-")) {  Proteins_DB <- Proteins_DB[-which(Proteins_DB$Sequence == "-"),] }
  if(any(Proteins_DB$Sequence == "No se encuentra el transcript id1" | Proteins_DB$Sequence == "No se encuentra el transcript id2")) {
    Proteins_DB <- Proteins_DB[-which(Proteins_DB$Sequence == "No se encuentra el transcript id1" | Proteins_DB$Sequence == "No se encuentra el transcript id2"),]
  }

  #Escribo el fasta con todas esas proteinas:
  fasta_file <- file(sprintf("%s/fasta_seq_proteins", path_dir), open = "w")
  writeLines(paste0(">", Proteins_DB$Transcript_id[1], "-", Proteins_DB$Gene[1]), fasta_file)
  # Escribir la secuencia correspondiente en la línea siguiente
  writeLines(Proteins_DB$Sequence[1], fasta_file)
  close(fasta_file)

  # Scan domains with HMMR 3
  pfam_path <- downloadPfam(omicsdo_sof_path = "/media/16TBDisk/Daniela/Pfam")
  system(sprintf("hmmscan --domtblout %s/found-domains.tab %s %s/fasta_seq_proteins", path_dir, pfam_path, path_dir))
  system(sprintf("cat %s/found-domains.tab | grep -v '^#' | sed 's/  */\t/g' | cut -f 1,2,4,20,21 > %s/found-domains-extract.tab", path_dir, path_dir))

  domains_df <- read.table(sprintf("%s/found-domains-extract.tab", path_dir), header = FALSE, sep = "\t",
                           col.names = c("Domain", "DomainID", "Gen", "Start", "End"))

  #domains_df <- domains_df[which(domains_df$Domain == "PK_Tyr_Ser-Thr"),]
  domains_df$TranscriptID <- sapply(strsplit(domains_df$Gen, split = "-"), `[`, 1)
  domains_df$Gen <- sapply(strsplit(domains_df$Gen, split = "-"), `[`, 2)

  domains_df$seq_peptide <- sequence_gene
  #domains_df$BP <- bp
  domains_df$Domain_Sequence <- substr(domains_df$seq_peptide, domains_df$Start, domains_df$End)


  #Conservacion de dominios: --------------------------------------
  d=1
  domains_df$conserved_domain <- "-"
  domains_df$porcentaje_conservacion <- "-"

  for (d in 1:nrow(domains_df)) {
    domain <- domains_df$Domain_Sequence[d]
    type_domain <- domains_df$Domain[d]
    message(sprintf("Analizando el porcentaje de conservación del dominio %s", type_domain))
    library(Biostrings)

    get_longest_common_substring <- function(domain, fusion) {
      domain_seq <- AAString(domain)
      fusion_seq <- AAString(fusion)
      hits <- pairwiseAlignment(domain_seq, fusion_seq, type = "local", scoreOnly = FALSE)
      conserved_domain <- as.character(subject(hits))
      return(conserved_domain)
    }
    conserved_domain <- get_longest_common_substring(domain = domain, fusion = seq_pedazo)
    domains_df$conserved_domain[d] <- conserved_domain

    porcentaje_conservacion <- round(nchar(conserved_domain) / nchar(domain) * 100)
    domains_df$porcentaje_conservacion[d] <- porcentaje_conservacion

    # Encontrar coordenadas del conserved_dominio dentro de la proteina:
    start_domain <- regexpr(conserved_domain, seq_prot)[1]
    end_domain <- start_domain + nchar(conserved_domain) -1
    domains_df$Coords_Conserved_Domain[d] <- paste0(start_domain, "-", end_domain)

    # Generar residuos para CBDOCK2:
    positions <- seq(start_domain, end_domain)
    amino_acids <- strsplit(domain, "")[[1]]
    domains_df$Residues_CBDOCK2[d] <- paste0(positions, ":", amino_acids, collapse = ",")
  }

  if(any(domains_df$Domain == "PK_Tyr_Ser-Thr")) {
    residuos <- domains_df$Residues_CBDOCK2[which(domains_df$Domain == "PK_Tyr_Ser-Thr")]
  } else {
    residuos <- "No Tirosin Kinase Domain"
  }

  domains_df <- domains_df[, c(1,10,2,3,4,5,6,7,8,9,11)]

  return(list(domains_df, residuos))
}

#-----------------------------------------------------------------------

# FUNCIONES PARAANOTADOR FUSIONES

source("/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/Funciones-AnotadorFusiones.R", echo=TRUE)

source("/media/16TBDisk/Daniela/OMICsdo/Grafico-Fusiones-Propio/draw_fusions_PROPIO.R", echo=TRUE)

#-----------------------------------------------------------------------

#BED_Anotador <- read_rds("/media/16TBDisk/Daniela/OMICsdo/data/BED_Anotador.RDS")

FLENI_ReporteFinal_QC <- read_excel("/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/FLENI_ReporteFinal_QC.xlsx")
FLENI_ReporteFinal_QC$Cohort <- "FLENI"
library(tidyr)
metrics_long <- FLENI_ReporteFinal_QC
metrics_long <- pivot_longer(
  metrics_long,
  cols = -c(ID, Cohort),
  names_to = "Métrica",
  values_to = "Valor"
)
metrics_long <- mutate(
  metrics_long,
  ID = as.character(ID),
  etiqueta = ID,
  color_id = ID
)
metrics_long_filtrado <- group_by(metrics_long, Métrica)
metrics_long_filtrado <- mutate(
  metrics_long_filtrado,
  Q1 = quantile(Valor, 0.25),
  Q3 = quantile(Valor, 0.75),
  IQR = Q3 - Q1,
  es_outlier = Valor < Q1 | Valor > Q3
)
unique(metrics_long_filtrado$Métrica)
metrics_long_filtrado$Fuente <- ifelse( metrics_long_filtrado$Métrica %in% c("Q_Mean_R1","%_>=Q30_R1","Q_Mean_R2","%_>=Q30_R2"), "QC Secuenciación", "QC Alineamiento")
colnames(metrics_long_filtrado)[3] <- "Metrica"
metrics_long_filtrado$Valor <- round(metrics_long_filtrado$Valor)

# SERVER
server <- function(input, output, session) {

  # Leer el Excel al iniciar la app
  datos <- reactive({
    read_excel(ruta_excel)
  })

  # DASHBOARD ----------------------------------------------------------------
  # Cantidad de muestras únicas
  output$num_muestras <- renderValueBox({
    n_muestras <- n_distinct(datos()$ID)  # suponiendo que la columna ID identifica muestras
    valueBox(
      value = n_muestras,
      subtitle = "Samples",
      icon = icon("user"),
      color = "blue"
    )
  })

  # Cantidad total de fusiones
  output$num_fusiones <- renderValueBox({
    n_fusiones <- nrow(datos())  # si cada fila es una fusión
    valueBox(
      value = n_fusiones,
      subtitle = "Fusions",
      icon = icon("dna"),
      color = "green"
    )
  })

  # Cantidad de cohortes únicas
  output$num_cohortes <- renderValueBox({
    n_cohortes <- n_distinct(datos()$Cohort)
    valueBox(
      value = n_cohortes,
      subtitle = "Cohorts",
      icon = icon("users"),
      color = "purple"
    )
  })

  # Cantidad de tipos de cáncer únicos
  output$num_tipos_cancer <- renderValueBox({
    n_cancer <- n_distinct(datos()$Cancer)
    valueBox(
      value = n_cancer,
      subtitle = "Types of Cancer",
      icon = icon("clinic-medical"),
      color = "red"
    )
  })

  # Gráfico de torta por cohorte (conteo de samples únicos)
  output$graficoCohorte <- renderPlotly({
    req(datos())

    resumen <- datos() %>%
      group_by(Cohort) %>%
      summarise(SamplesUnicos = n_distinct(ID), .groups = "drop")

    plot_ly(
      data = resumen,
      labels = ~Cohort,
      values = ~SamplesUnicos,
      type = "pie",
      textinfo = "label+percent",
      hoverinfo = "label+value+percent"
    ) %>%
      layout(title = "Samples per Cohort")
  })

  # Gráfico de torta por cancer (conteo de samples únicos)
  output$graficoCancer <- renderPlotly({
    req(datos())

    resumen <- datos() %>%
      group_by(Cancer) %>%
      summarise(SamplesUnicos = n_distinct(ID), .groups = "drop")

    plot_ly(
      data = resumen,
      labels = ~Cancer,
      values = ~SamplesUnicos,
      type = "pie",
      textinfo = "label+percent",
      hoverinfo = "label+value+percent"
    ) %>%
      layout(title = "Samples per Cancer type")
  })

  # --- Gráfico de barras: cantidad de fusiones por Cohorte ---
  output$FusCohorte <- renderPlotly({
    req(datos())

    resumen_cohorte <- datos() %>%
      group_by(Cohort) %>%
      summarise(Fusiones = n(), .groups = "drop")

    plot_ly(
      data = resumen_cohorte,
      x = ~Cohort,
      y = ~Fusiones,
      type = "bar",
      text = ~Fusiones,
      textposition = "auto",
      hoverinfo = "x+y"
    ) %>%
      layout(
        title = "Nº of Fusions per Cohort",
        xaxis = list(title = "Cohort"),
        yaxis = list(title = "Nº of Fusions")
      )
  })


  # --- Gráfico de barras: cantidad de fusiones por Cancer ---
  output$FusCancer <- renderPlotly({
    req(datos())

    resumen_cancer <- datos() %>%
      group_by(Cancer) %>%
      summarise(Fusiones = n(), .groups = "drop")

    plot_ly(
      data = resumen_cancer,
      x = ~Cancer,
      y = ~Fusiones,
      type = "bar",
      text = ~Fusiones,
      textposition = "auto",
      hoverinfo = "x+y"
    ) %>%
      layout(
        title = "Nº Fusions per Cancer Type",
        xaxis = list(title = "Cancer"),
        yaxis = list(title = "Nº of Fusions")
      )
  })

  # Boxplot 1: Número de fusiones por ID, agrupado por Cancer

  output$seleccionGen <- renderUI({
    req(input$columnaFiltrarGen)
    df <- datos()
    genes_posibles <- unique(df[[input$columnaFiltrarGen]])

    selectInput("genesSeleccionados", "Selecciona genes:",
                choices = genes_posibles,
                multiple = TRUE,
                selected = genes_posibles[1])
  })

  datosAgrupados <- reactive({
    req(input$columnasAgrupar)
    df <- datos()
    # Filtramos primero por genes seleccionados
    if(!is.null(input$genesSeleccionados)){
      df <- df %>% filter(.data[[input$columnaFiltrarGen]] %in% input$genesSeleccionados)
    }
    # Luego agrupamos y resumimos
    df %>%
      group_by(ID, across(all_of(input$columnasAgrupar))) %>%
      summarise(Fusiones = n(), .groups = "drop")
  })

  output$BoxplotCancer <- renderPlotly({
    df <- datos()
    req(df)
    req(input$columnasAgrupar)

    columnas_label <- setdiff(input$columnasAgrupar, "Cohort")

    df$contig1 <- factor(df$contig1, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
                                                   "17", "18", "19", "20", "X", "Y"))
    df$contig2 <- factor(df$contig2, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
                                                "17", "18", "19", "20", "X", "Y"))
    df$confidence <- factor(df$confidence, levels = c("low", "medium", "high"))

    # Filtrado por genes seleccionados (si hay)
    if(!is.null(input$genesSeleccionados) && length(input$genesSeleccionados) > 0){
      df <- df %>% filter(.data[[input$columnaFiltrarGen]] %in% input$genesSeleccionados)
      incluir_gen <- TRUE
    } else {
      incluir_gen <- FALSE
    }

    # Creamos la etiqueta para el eje X
    if(length(columnas_label) > 0){
      df <- df %>%
        mutate(Xaxis = apply(across(all_of(columnas_label)), 1, paste, collapse = " - "))

      # Crear los niveles ordenados para Xaxis
      niveles_ordenados <- expand.grid(
        lapply(columnas_label, function(col) {
          if(col %in% c("contig1", "contig2")) {
            levels(df[[col]])[levels(df[[col]]) %in% unique(df[[col]])]
          } else if(col == "confidence") {
            levels(df[[col]])[levels(df[[col]]) %in% unique(df[[col]])]
          } else {
            sort(unique(df[[col]]))
          }
        })
      )

      # Invertir el orden de las columnas para expand.grid
      nombres_cols <- rev(columnas_label)
      niveles_ordenados <- expand.grid(
        setNames(
          lapply(rev(columnas_label), function(col) {
            if(col %in% c("contig1", "contig2")) {
              levels(df[[col]])[levels(df[[col]]) %in% unique(df[[col]])]
            } else if(col == "confidence") {
              levels(df[[col]])[levels(df[[col]]) %in% unique(df[[col]])]
            } else {
              sort(unique(df[[col]]))
            }
          }),
          rev(columnas_label)
        )
      )

      # Crear la cadena de niveles en el mismo orden que Xaxis
      niveles_xaxis <- apply(niveles_ordenados[, rev(names(niveles_ordenados)), drop = FALSE], 1, paste, collapse = " - ")
      niveles_existentes <- niveles_xaxis[niveles_xaxis %in% unique(df$Xaxis)]
      df$Xaxis <- factor(df$Xaxis, levels = niveles_existentes)


    } else {
      df <- df %>% mutate(Xaxis = "")
    }

    # Si hay genes seleccionados, los agregamos al eje X
    if(incluir_gen){
      df <- df %>% mutate(Xaxis = paste(Xaxis, .data[[input$columnaFiltrarGen]], sep = " | "))
      # Reordenar considerando también los genes
      df$Xaxis <- factor(df$Xaxis, levels = sort(unique(df$Xaxis)))
    }

    # Convertimos la columna de color a character
    color_var <- if("Cohort" %in% input$columnasAgrupar) "Cohort" else input$columnasAgrupar[1]
    df[[color_var]] <- as.character(df[[color_var]])


    if(isTRUE(input$usarTotalFusiones)){
      # Gráfico de barras con total de fusiones
      df_total <- df %>%
        group_by(Xaxis, .data[[color_var]]) %>%
        summarise(TotalFusiones = n(), .groups = "drop")


      plot_ly(df_total,
              x = ~Xaxis,
              y = ~TotalFusiones,
              type = "bar",
              color = ~.data[[color_var]],
              colors = "Set2") %>%
        layout(
          xaxis = list(title = if(length(columnas_label) > 0) paste(columnas_label, collapse = " + ") else "Cohort",
                       categoryorder = "array",
                       categoryarray = levels(df_total$Xaxis)),
          yaxis = list(title = "Nº of fusions"),
          barmode = "group",
          legend = list(title = list(text = color_var))
        )

    } else {
      # ------------- Boxplot por ID ---------------------------
      df_box <- df %>%
        group_by(ID, Xaxis, .data[[color_var]]) %>%
        summarise(Fusiones = n(), .groups = "drop")

      plot_ly(df_box,
              x = ~Xaxis,
              y = ~Fusiones,
              type = "box",
              color = ~.data[[color_var]],
              colors = "Set2",
              boxpoints = "all",
              jitter = 0.5,
              pointpos = 0) %>%
        layout(
          xaxis = list(
            title = if(length(columnas_label) > 0) paste(columnas_label, collapse = " + ") else "Cohort",
            showticklabels = TRUE,
            categoryorder = "array",
            categoryarray = levels(df_box$Xaxis)
          ),
          yaxis = list(title = "Nº of Fusions per Sample"),
          boxmode = "group",
          legend = list(title = list(text = color_var))
        )
    }
  })


  #Barras con TK:
  output$TK_100_Cohorte <- renderPlotly({
    req(datos())

    df_counts <- datos() %>%
      mutate(Kinase = ifelse(grepl("Protein_kinase_domain\\(100%\\)", retained_protein_domains),
                             "Kinase domain", "Other")) %>%
      group_by(Cohort, Kinase) %>%
      summarise(Fusions = n(), .groups = "drop")

    p <- ggplot(df_counts, aes(x = Cohort, y = Fusions, fill = Kinase)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = "100% retained Kinase domain",
           x = "Cohort", y = "Fusion Count") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggplotly(p)
  })

  output$TK_Cohorte <- renderPlotly({
    req(datos())

    df_counts <- datos() %>%
      mutate(Kinase = ifelse(grepl("Protein_kinase_domain", retained_protein_domains),
                             "Kinase domain", "Other")) %>%
      group_by(Cohort, Kinase) %>%
      summarise(Fusions = n(), .groups = "drop")

    p <- ggplot(df_counts, aes(x = Cohort, y = Fusions, fill = Kinase)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = "Any % Retained Kinase Domain",
           x = "Cohort", y = "Fusion Count") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggplotly(p)
  })

  # Mostrar tabla
  output$tabla <- renderDT({
    datatable(datos(), options = list(pageLength = 10, scrollX = TRUE))
  })


  # BY PATIENT ---------------------------

  observe({
    updateSelectInput(session, "selectCohorte",
                      choices = sort(unique(datos()$Cohort)))
  })

  # Actualizar selectInput ID en función del Cohorte elegido
  observeEvent(input$selectCohorte, {
    ids <- datos() %>% filter(Cohort == input$selectCohorte) %>% pull(ID) %>% unique()
    updateSelectInput(session, "selectID",
                      choices = c("All" = "", sort(ids)),
                      selected = "")
  })
  # Filtrar tabla según Cohorte, ID y fusiones conocidas
  datosPaciente <- reactive({
    req(input$selectCohorte)

    if(!is.null(input$selectID) && input$selectID != "") {
      df <- datos() %>%
      filter(Cohort == input$selectCohorte,
             ID == input$selectID)
    } else {
      df <- datos() %>%
        filter(Cohort == input$selectCohorte)
    }

    # Filtrar solo fusiones conocidas si la casilla está tildada
    if(!is.null(input$fusionesConocidas) && input$fusionesConocidas){
      genesConocidos <- c("BRAF", "NTRK1", "NTRK2", "NTRK3", "RET", "MAPK", "ROS1",
                          "ALK", "BCOR", "CIC", "EGFR", "ETV1", "EWSR1", "FGFR1", "FGFR2", "FGFR3", "FOXR2", "FUS",
                          "HTRA1", "MET", "MGMT", "MN1", "MYB", "MYBL1", "MYC", "PCSK5", "PDGFRA", "PKD1", "PVT1",
                          "RAF1", "RECK", "YAP1", "ZFTA(C11orf95)", "ERBB2")
      df <- df %>%
        filter(gene1 %in% genesConocidos | gene2 %in% genesConocidos)
    }

    # Filtrar solo fusiones tirosin kinasa
    if(!is.null(input$fusionesTK) && input$fusionesTK){
      df <- df %>%
        filter(grepl("Protein_kinase_domain", retained_protein_domains))
    }

    return(df)
  })


  # Renderizar tabla manteniendo selección múltiple y scroll
  output$tablaPaciente <- renderDT({
    datatable(
      datosPaciente(),
      selection = "multiple",
      options = list(scrollX = TRUE)
    )
  })


  # Obtener la fila seleccionada
  fusionesSeleccionadas <- reactive({
    req(input$tablaPaciente_rows_selected)
    datosPaciente()[input$tablaPaciente_rows_selected, ]
  })

  ###########################
  output$fusionPanels <- renderUI({
    req(fusionesSeleccionadas())

    tabs <- lapply(1:nrow(fusionesSeleccionadas()), function(i) {
      fusion <- fusionesSeleccionadas()[i, ]
      fusion_id <- paste0("fusion_", i)  # id único

      tabPanel(
        title = paste("Fusion:", fusion$gene1, "-", fusion$gene2),
        tabsetPanel(
          tabPanel("Fusion plot", plotOutput(paste0("fusionPlot_", fusion_id))),
          tabPanel("Domains",
                   # Gráfico arriba ocupando todo el ancho
                   fluidRow(column(12, plotOutput(paste0("domains_plot_", fusion_id))) ),
                   # Tablas lado a lado
                   fluidRow( column(6, DTOutput(paste0("domains_gene1_", fusion_id))),
                     column(6,  DTOutput(paste0("domains_gene2_", fusion_id))))
                   ),
          tabPanel("Protein sequence", uiOutput(paste0("fusionSeq_", fusion_id))),
          tabPanel("For docking",
                   tabsetPanel(
                     tabPanel("Gene1", verbatimTextOutput(paste0("fusionDockingGene1_", fusion_id))),
                     tabPanel("Gene2", verbatimTextOutput(paste0("fusionDockingGene2_", fusion_id)))
                   )
          ),
          tabPanel("Drugs",
                   tabsetPanel(
                     tabPanel("Gene 1 Drugs", DTOutput(paste0("drugs_gene1_", fusion_id))),
                     tabPanel("Gene 2 Drugs", DTOutput(paste0("drugs_gene2_", fusion_id)))
                   )
          )
        )
      )
    })

    do.call(tabBox, c(width = 12, tabs))
  })

  observe({
    req(fusionesSeleccionadas())

    for (i in 1:nrow(fusionesSeleccionadas())) {
      fusion <- fusionesSeleccionadas()[i, ]
      fusion_id <- paste0("fusion_", i)  # sufijo único

      # --- Subventana 1: gráfico de la fusión ---
      local({
        output[[paste0("fusionPlot_", fusion_id)]] <- renderPlot({
          f <- fusion; id <- fusion_id
          # Validaciones visibles en la UI (si falta algo, lo muestra)
          validate(
            need(!is.null(f), "fusion is NULL"),
            need(!grepl(",", f$gene1), "2 genes in position 5'"),
            need(!grepl(",", f$gene2), "2 genes in position 3'"),
            need((f$peptide_sequence != "."), "This fusion does not generate a protein."),
            need(exists("exons"), "exons object missing"),
            need(exists("plot_Fusion_without_coverage"), "plot_Fusion_without_coverage() not found")
          )

          # Caso sin proteína: dibuja mensaje en vez de dejar vacío
          #if( f$peptide_sequence == ".") {
          #  plot.new()
          #  text(0.5, 0.5, "This fusion cannot be transcribed into a protein", cex = 1.2)
          #  return()
          #}

          print(plot_Fusion_without_coverage(fusions = f, exons = exons))


          #plot_Fusion(fusions = datos[1,])
        })
      })


      # --- Subventana 2: dominios fusion: 2 tablas y 1 gráfico---
      local({
        f <- fusion; id <- fusion_id

        if(f$peptide_sequence != "." & !grepl(",", f$gene1) & !grepl(",", f$gene2)) {
          seq_total <- AnotateFusion(
                transcript_id1   = f$transcript_id1,
                transcript_id2   = f$transcript_id2,
                protein_df       = protein_df,
                peptide_sequence = f$peptide_sequence
              )
              seq_prot      <- seq_total[[1]]
              sequence_gene1 <- seq_total[[2]]
              sequence_gene2 <- seq_total[[3]]
              seq_pedazo_1  <- seq_total[[4]]
              seq_pedazo_2  <- seq_total[[5]]

              if(seq_prot == ".NO SE ENCUENTRA EL PEDAZO DE PEPTIDE EN LA SECUENCIA DEL GEN") {
                seq_prot      <- "This fusion cannot be transcribed into a protein"
                sequence_gene1 <- "-"
                sequence_gene2 <- "-"
                seq_pedazo_1  <- "-"
                seq_pedazo_2  <- "-"

                domains_df1 <- data.frame( Message =  "This fusion cannot be transcribed into a protein")
                domains_df2 <- data.frame( Message =  "This fusion cannot be transcribed into a protein")

              } else {
                domains_df1 <- GetDomain(
                  transcript_id = f$transcript_id1,
                  gene          = f$gene1,
                  sequence_gene = sequence_gene1,
                  seq_pedazo    = seq_pedazo_1,
                  seq_prot      = seq_prot
                )
                domains_df1 <- domains_df1[[1]]

                domains_df2 <- GetDomain(
                  transcript_id = f$transcript_id2,
                  gene          = f$gene2,
                  sequence_gene = sequence_gene2,
                  seq_pedazo    = seq_pedazo_2,
                  seq_prot      = seq_prot
                )
                domains_df2 <- domains_df2[[1]]

              }

        } else {
          domains_df1 <- data.frame(Message = "This fusion does not generate a protein.")
          domains_df2 <- data.frame(Message = "This fusion does not generate a protein.")
        }

        # --- Tablas de dominios ---
        output[[paste0("domains_gene1_", id)]] <- DT::renderDataTable({
          domains_df1
        }, options = list(scrollX = TRUE, scrollY = "400px", paging = FALSE))

        output[[paste0("domains_gene2_", id)]] <- DT::renderDataTable({
          domains_df2
        }, options = list(scrollX = TRUE, scrollY = "400px", paging = FALSE))

        # --- Gráfico de dominios ---

        output[[paste0("domains_plot_", id)]] <- renderPlot({
          f <- fusion; id <- fusion_id
          # No hacer nada si no hay proteína
          if(f$peptide_sequence == ".") return(NULL)

          # Paso 1: Preparar datos para ggplot
          library(ggplot2)
          library(dplyr)
          library(stringr)

          # Función para preparar los dominios y apilarlos verticalmente
          prepare_domains_plot <- function(domains_df, y_bottom, height=0.2) {
            n_domains <- nrow(domains_df)
            if(n_domains == 0) return(NULL)

            # Parsear Coords_Conserved_Domain tipo "20-22" y ordenar
            domains_df <- domains_df %>%
              mutate(
                start_coord = as.numeric(str_extract(Coords_Conserved_Domain, "^[0-9]+")),
                end_coord   = as.numeric(str_extract(Coords_Conserved_Domain, "(?<=-)[0-9]+"))
              ) %>%
              #arrange(start_coord) %>%  # <- Ordenar por posición inicial
              mutate(
                ymin = y_bottom + seq(0, by=height, length.out = n_domains),
                ymax = ymin + height
              )
            return(domains_df)
          }

          gene_length1 <- nchar(seq_pedazo_1)
          gene_length2 <- nchar(seq_pedazo_2)

          # Preparar dominios para ggplot
          domains_plot1 <- prepare_domains_plot(domains_df1, y_bottom=1)
          domains_plot2 <- prepare_domains_plot(domains_df2, y_bottom=1)

          domains_plot1 <- domains_plot1 %>% mutate(xmin = start_coord, xmax = end_coord)
          domains_plot2 <- domains_plot2 %>% mutate(xmin = start_coord, xmax = end_coord)

          max_domains <- max(nrow(domains_plot1), nrow(domains_plot2), 1)
          gene_height <- max_domains * 0.2 + 0.1

          gene_rects <- data.frame(
            gene = c(domains_df1$Gen[1], domains_df2$Gen[1]),
            xmin = c(0, nchar(seq_pedazo_1)),
            xmax = c(nchar(seq_pedazo_1), nchar(seq_pedazo_1) + nchar(seq_pedazo_2)),
            ymin = c(0.9, 0.9),
            ymax = c(0.9 + gene_height, 0.9 + gene_height),
            fill_color = c("#DCD0FF", "#FFDBBB")
          )

          ggplot() +
            geom_rect(data=gene_rects,
                      aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                      fill=gene_rects$fill_color) +
            geom_rect(data=domains_plot1,
                      aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Domain)) +
            geom_rect(data=domains_plot2,
                      aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=Domain)) +
            geom_text(data=domains_plot1,
                      aes(x=(xmin+xmax)/2, y=(ymin+ymax)/2,
                          label=paste0(porcentaje_conservacion,"%")),
                      size=3, color="black") +
            geom_text(data=domains_plot2,
                      aes(x=(xmin+xmax)/2, y=(ymin+ymax)/2,
                          label=paste0(porcentaje_conservacion,"%")),
                      size=3, color="black") +
            geom_text(data=gene_rects,
                      aes(x=(xmin+xmax)/2, y=ymax + 0.05,
                          label=paste0(gene, " (", xmax - xmin, " AA)")),
                      size=4, fontface="bold", vjust=0) +
            scale_y_continuous(breaks = NULL) +
            theme_minimal() +
            xlab("Amino acid position") + ylab("") +
            theme(legend.position="left")
        })
      })


      # --- Subventana 3: secuencia proteína ---
      local({
        f <- fusion; id <- fusion_id
        output[[paste0("fusionSeq_", id)]] <- renderUI({

          seq_total <- AnotateFusion(
            transcript_id1   = f$transcript_id1,
            transcript_id2   = f$transcript_id2,
            protein_df       = protein_df,
            peptide_sequence = f$peptide_sequence
          )

          seq_prot <- seq_total[[1]]
          sequence1 <- seq_total[[2]]
          sequence2 <- seq_total[[3]]

          if( grepl(",", f$gene1) | grepl(",", f$gene2) |  f$peptide_sequence == "." | seq_prot == ".NO SE ENCUENTRA EL PEDAZO DE PEPTIDE EN LA SECUENCIA DEL GEN") {
            seq_prot <- "This fusion cannot be transcribed into a protein"
            sequence1 <- "-"
            sequence2 <- "-"
          }

          HTML(
            paste0(
              "<h4>Protein sequence resulted from the fusion:</h4>",
              "<pre>", seq_prot, "</pre>",
              "<h4>Gene 1 Complete Peptide Sequence:</h4>",
              "<pre>", sequence1, "</pre>",
              "<h4>Gene 2 Complete Peptide Sequence:</h4>",
              "<pre>", sequence2, "</pre>"
            )
          )
        })
      })


      # --- Subventana 4: docking Gene1 ---
      local({
        f <- fusion; id <- fusion_id
        output[[paste0("fusionDockingGene1_", id)]] <- renderPrint({
          seq_total <- AnotateFusion(
            transcript_id1   = f$transcript_id1,
            transcript_id2   = f$transcript_id2,
            protein_df       = protein_df,
            peptide_sequence = f$peptide_sequence
          )
          seq_prot <- seq_total[[1]]
          sequence_gene1 <- seq_total[[2]]
          seq_pedazo_1<- seq_total[[4]]
          d <- GetDomain(
            transcript_id = f$transcript_id1,
            gene          = f$gene1,
            sequence_gene = sequence_gene1,
            seq_pedazo = seq_pedazo_1,
            seq_prot      = seq_prot
          )
          d[[2]]  # residuos
        })
      })

      # --- Subventana 4: docking Gene2 ---
      local({
        f <- fusion; id <- fusion_id
        output[[paste0("fusionDockingGene2_", id)]] <- renderPrint({
          seq_total <- AnotateFusion(
            transcript_id1   = f$transcript_id1,
            transcript_id2   = f$transcript_id2,
            protein_df       = protein_df,
            peptide_sequence = f$peptide_sequence
          )
          seq_prot <- seq_total[[1]]
          sequence_gene2 <- seq_total[[3]]
          seq_pedazo_2 <- seq_total[[5]]
          d <- GetDomain(
            transcript_id = f$transcript_id2,
            gene          = f$gene2,
            sequence_gene = sequence_gene2,
            seq_pedazo = seq_pedazo_2,
            seq_prot      = seq_prot
          )
          d[[2]]  # residuos
        })
      })

      # --- Subventana 5: drogas gene1 ---
      local({
        f <- fusion; id <- fusion_id
        output[[paste0("drugs_gene1_", id)]] <- DT::renderDataTable({
          getDrugs(protein_df = protein_df, gene = f$gene1)
        }, options = list(scrollX = TRUE, scrollY = "400px", paging = FALSE))
      })

      # --- Subventana 5: drogas gene2 ---
      local({
        f <- fusion; id <- fusion_id
        output[[paste0("drugs_gene2_", id)]] <- DT::renderDataTable({
          getDrugs(protein_df = protein_df, gene = f$gene2)
        }, options = list(scrollX = TRUE, scrollY = "400px", paging = FALSE))
      })
    }
  })

  # Metadata by Patient ------------------
  datos_meta <- reactive({
    Metadata_FLENI
  })

  datosPacienteMeta <- reactive({
    req(input$selectCohorte, input$selectID)

    df_meta <- datos_meta() %>%
      filter(Cohort == input$selectCohorte,
             ID == input$selectID)
    return(df_meta)
  })
  # Renderizar tabla
  output$tablaPacienteMeta <- renderDT({
    datatable(
      #datosPacienteMeta()[,-which(is.na(datosPacienteMeta()))],  # Solo columnas que querés mostrar
      datosPacienteMeta(),
      options = list(scrollX = TRUE)
    )
  })

  # QC by Patient ------------------
  datos_qc <- reactive({
    metrics_long_filtrado
  })

  datosPacienteQC <- reactive({
    req(input$selectCohorte, input$selectID)

    df <- datos_qc() %>%
      filter(Cohort == input$selectCohorte,
             ID == input$selectID)
    return(df)
  })
  # Renderizar tabla
  output$tablaPacienteQC <- renderDT({
    datatable(
      datosPacienteQC()[, c("Metrica", "Valor")],  # Solo columnas que querés mostrar
      options = list(scrollX = TRUE)
    )
  })


  output$QCBarras <- renderPlotly({
    req(input$selectCohorte, input$selectID)

    df <- datos_qc() %>%
      filter(Cohort == input$selectCohorte)

    # Calcular percentiles
    #df_percentil <- df %>%
    #  group_by(Metrica) %>%
    #  mutate(Percentil = percent_rank(Valor) * 100) %>%
    #  ungroup()

    #score de centralidad en el boxplot:
    df_score <- df %>%
      group_by(Metrica) %>%
      mutate(
        Q1 = quantile(Valor, 0.25, na.rm = TRUE),
        Q3 = quantile(Valor, 0.75, na.rm = TRUE),
        Mediana = median(Valor, na.rm = TRUE),
        IQR = Q3 - Q1,
        Score = case_when(
          # Todos los valores iguales
          IQR == 0 ~ 100,
          # Dentro del box
          Valor >= Q1 & Valor <= Q3 ~ 70 + 30 * (1 - abs(Valor - Mediana) / (IQR / 2)),
          # Moderadamente por debajo (entre Q1 - IQR y Q1)
          Valor < Q1 & Valor >= (Q1 - IQR) ~ 30 + 40 * (Valor - (Q1 - IQR)) / IQR,
          # Muy por debajo (menor a Q1 - 5*IQR)
          Valor < (Q1 - 5*IQR) ~ pmax(0, 30 * (1 - (Q1 - Valor) / (10 * IQR))),
          # Moderadamente por debajo (entre Q1 - 5*IQR y Q1 - IQR)
          Valor < (Q1 - IQR) & Valor >= (Q1 - 5*IQR) ~ 30 + 40 * (Valor - (Q1 - 5*IQR)) / (4 * IQR),
          # Moderadamente por encima (entre Q3 y Q3 + 5*IQR)
          Valor > Q3 & Valor <= (Q3 + 5*IQR) ~ 30 + 40 * ((Q3 + 5*IQR) - Valor) / (5 * IQR),
          # Muy por encima (mayor a Q3 + 5*IQR)
          Valor > (Q3 + 5*IQR) ~ pmax(0, 30 * (1 - (Valor - Q3) / (10 * IQR))),
          # Fallback
          TRUE ~ 0
        )
      ) %>%
      ungroup()


    # Datos del ID seleccionado
    df_id <- df_score[which(df_score$ID == input$selectID),]
    metricas <- unique(df_score$Metrica)

    # Gráfico vacío
    p <- plot_ly()
    list_plots <- list()

    for (i in 1:length(metricas)) {
      m <- metricas[i]
      # Valor de score del ID
      x_val <- df_id$Score[i]

      p <- plot_ly() %>%
        layout(
          shapes = list(
            list(type = "rect", x0 = 0, x1 = 30, y0 = 0, y1 = 1,
                 fillcolor = "red", line = list(width = 0)),
            list(type = "rect", x0 = 30, x1 = 70, y0 = 0, y1 = 1,
                 fillcolor = "yellow", line = list(width = 0)),
            list(type = "rect", x0 = 70, x1 = 100, y0 = 0, y1 = 1,
                 fillcolor = "green", line = list(width = 0)),
            list(type = "rect", x0 = x_val-0.2, x1 = x_val+0.2, y0 = 0, y1 = 1,
                 fillcolor = "black", line = list(width = 0))
          ),
          xaxis = list(range = c(-5, 110), title = "Score"),
          yaxis = list(range = c(0, 1.5),
                       showticklabels = FALSE,
                       showgrid = FALSE,
                       zeroline = FALSE),
          margin = list(l = 50, r = 50, t = 30, b = 30)
        ) %>%
        add_trace(
          x = x_val,
          y = 1.5,  # 🔹 encima de la barra
          type = "scatter",
          mode = "markers",
          marker = list(symbol = "triangle-down", size = 20, color = "black"),
          text = paste0("Métrica: ", m,
                        "<br>ID: ", unique(df_id$ID),
                        "<br>Score: ", round(x_val, 1)),
          hoverinfo = "text",
          showlegend = FALSE
        )

      list_plots[[i]] <- p

    }
    # Combinar todos los plots en subplots verticales
    subplot(list_plots, nrows = length(list_plots), shareX = TRUE, titleY = TRUE)
  })


  output$explicacion_score_qc <- renderUI({
    HTML("<b>Interpretación del score:</b><br>
       El score refleja cuán cercano está ese ID a los valores de las demás muestras.<br>
       <ul>
         <li><b>70–100:</b> valor dentro del rango esperado (entre Q1 y Q3), siendo 100 la mediana.</li>
         <li><b>30–70:</b> valor fuera del rango intercuartílico, pero todavía cercano (moderado desvío).</li>
         <li><b>0–30:</b> valor extremadamente alejado (muy por debajo de Q1−5·IQR o muy por encima de Q3+5·IQR).</li>
       </ul>")
  })

  # REPORT: -------------------------------------

  output$downloadReport <- downloadHandler(
    filename = function() {
      paste0("FusionsReport_", input$selectID, ".docx")
    },
    content = function(file) {
      # Paquetes para escribir DOCX
      library(officer)
      library(flextable)

      doc <- read_docx()

      # 🔹 Definir el encabezado como un block_list
      header_content <- block_list(
        fpar(ftext(paste("Nombre del paciente:", input$selectID),
                   prop = fp_text(font.size = 12, bold = TRUE))),
        fpar(ftext(paste("Número de protocolo:", input$selectID),
                   prop = fp_text(font.size = 12, bold = TRUE))),
        fpar(ftext(paste("Fecha de entrega del informe:",
                         format(Sys.Date(), "%d/%m/%Y")),
                   prop = fp_text(font.size = 10)))
      )

      # 🔹 Reemplazar encabezado en todas las secciones
      # CORREGIDO: Usar body_add_block_list() para agregar el objeto header_content
      doc <- doc %>%
        body_add_fpar(header_content[[1]]) %>%
        body_add_fpar(header_content[[2]]) %>%
        body_add_fpar(header_content[[3]])


      # 🔹 Definir la tabla
      info_df <- data.frame(
        "Número de Protocolo" = "36919",
        "Nombre" = "---",
        "Fecha de recepción de la muestra" = "01/07/2019",
        "Tipo de muestra" = "---",
        stringsAsFactors = FALSE
      )

      # 🔹 Convertir en flextable
      info_ft <- flextable(info_df)

      # 🔹 Agregar al documento
      doc <- doc %>%
        body_add_par("PACIENTE", style = "heading 2") %>%
        body_add_flextable(info_ft) %>%
        body_add_par(" ", style = "Normal")


      doc <- doc %>%
        body_add_par("ANÁLISIS", style = "heading 2") %>%
        body_add_par("Tipo de estudio", style = "heading 3") %>%
        body_add_par("Secuenciación masiva de ARN total obtenido de tejido tumoral congelado.", style = "Normal") %>%
        body_add_par("Metodología", style = "heading 3") %>%
        body_add_par("Las muestras se prepararon de acuerdo a las guías de preparación del kit TruSeq Stranded Total RNA LT. Las bibliotecas fueron secuenciadas con el secuenciador de la plataforma NovaSeq 6000.", style = "Normal") %>%
        body_add_par("Softwares utilizados", style = "heading 3") %>%
        body_add_par(paste(
          "Para el análisis de los datos de secuenciación se utilizaron los siguientes programas:",
          "1. FastQC v0.11.7 (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)",
          "2. TrimGalore v0.6.10",
          "3. STAR v2.7.11a",
          "4. ARRIBA v2.4.0",
          "Para realizar el mapeo de las lecturas se utilizó como referencia el hg38 de Ensembl.",
          sep = "\n"),
          style = "Normal")

      # QC info
      qc_df <- datosPacienteQC()
      qc_df <- as.data.frame(qc_df[, c(3,4,11)])
      doc <- doc %>%
        body_add_par("Quality Metrics", style = "heading 2") %>%
        body_add_flextable(flextable(qc_df)) %>%
        body_add_par(" ", style = "Normal")

      # 🔹 Fusiones seleccionadas
      fusiones <- fusionesSeleccionadas()
      doc <- doc %>%
        body_add_par("RESULTADOS", style = "heading 3")

      # CORREGIDO: la condición del if para evitar un error lógico
      if(nrow(fusiones) > 0) {
        for (i in 1:nrow(fusiones)) {
          f <- fusiones[i, ]
          doc <- doc %>%
            body_add_par(paste0("Fusion ", i, ": ", f$gene1, " - ", f$gene2), style = "heading 2") %>%
            body_add_par(paste("Transcript 1:", f$transcript_id1), style = "Normal") %>%
            body_add_par(paste("Transcript 2:", f$transcript_id2), style = "Normal") %>%
            body_add_par(paste("Peptide sequence:", f$peptide_sequence), style = "Normal")

          # Opcional: agregar tabla de drogas asociadas
          drugs1 <- getDrugs(protein_df = protein_df, gene = f$gene1)
          drugs2 <- getDrugs(protein_df = protein_df, gene = f$gene2)

          doc <- doc %>%
            body_add_par("Drugs for Gene 1", style = "heading 3") %>%
            body_add_flextable(flextable(drugs1)) %>%
            body_add_par("Drugs for Gene 2", style = "heading 3") %>%
            body_add_flextable(flextable(drugs2)) %>%
            body_add_par(" ", style = "Normal")
        }
      } else {
        doc <- doc %>%
          body_add_par("El análisis de fusiones mediante secuenciación de ARN no reveló fusiones en los genes de interés.", style = "Normal")
      }

      doc <- doc %>%
        body_add_par("Los datos de secuenciación quedarán disponibles para un posible pedido de re-análisis y búsqueda de fusiones en otros genes.", style = "Normal")

      doc <- doc %>%
        body_add_par("INFORME", style = "heading 3") %>%
        body_add_par("Interpretación Clínica", style = "heading 2") %>%
        body_add_par(" ", style = "Normal") %>%
        body_add_par("Limitaciones de la técnica", style = "heading 2") %>%
        body_add_par("El análisis de fusiones mediante secuenciación de ARN no detecta inserciones o deleciones grandes, ni cambios en el número de copias (CNVs). Un resultado negativo no descarta la presencia de fusiones o rearreglos que estén por debajo del límite de detección de la prueba. La presencia de polimorfismos atípicos puede llevar a resultados que sean falsos positivos o negativos.", style = "Normal") %>%
        body_add_par(" ", style = "Normal") %>%
        body_add_par("Información adicional", style = "heading 2") %>%
        body_add_par(paste0("Ensayos clínicos disponibles al día de la fecha pueden ser encontrados en los siguientes sitios:",
                            "1. ClinicalTrials.gov (www.clinicaltrials.gov)",
                            "2. Mayo Clinic (www.mayo.edu/research/clinical-trials)",
                            "3. National Cancer Institute (www.cancer.gov/clinicaltrials)", sep = "\n"), style = "Normal") %>%
        body_add_par("Referencias", style = "heading 2")

      # Guardar el documento
      print(doc, target = file)
    }
  )

  # SERVER: by Gene ----------------------------------------------
  datos_genes <- reactive({
    read_excel(ruta_excel_genes)
  })

  observe({
    updateSelectInput(session, "selectCohorteGenes",  # ID único para la sección genes
                      choices = sort(unique(datos_genes()$Cohort)))
  })

  # Filtrar tabla según Cohorte y fusiones conocidas
  datosGenes <- reactive({
    req(input$selectCohorteGenes)  # Usar el ID correspondiente
    df <- datos_genes() %>% filter(Cohort == input$selectCohorteGenes)

    # Filtrar solo fusiones conocidas si la casilla está tildada
    if(!is.null(input$genesConocidos) && input$genesConocidos){
      genesC <- c("BRAF", "NTRK1", "NTRK2", "NTRK3", "RET", "MAPK", "ROS1",
                          "ALK", "BCOR", "CIC", "EGFR", "ETV1", "EWSR1", "FGFR1", "FGFR2", "FGFR3", "FOXR2", "FUS",
                          "HTRA1", "MET", "MGMT", "MN1", "MYB", "MYBL1", "MYC", "PCSK5", "PDGFRA", "PKD1", "PVT1",
                          "RAF1", "RECK", "YAP1", "ZFTA(C11orf95)", "ERBB2")
      df <- df %>%
        filter(Gene %in% genesC)
    }
    return(df)
  })

  # Renderizar la tabla de genes
  output$tablaGenes <- renderDT({
    datatable(datosGenes(),
              options = list(scrollX = TRUE, pageLength = 10),
              filter = 'top')
  })
  # SERVER: by Fusion --------------------------------------------

  observeEvent(input$gen5, {
    req(input$gen5)

    bp_choices <- unique(datos()$breakpoint1[which(datos()$gene1 == input$gen5)])

    # Si no hay coincidencias → usar vector vacío
    if (length(bp_choices) == 0 || all(is.na(bp_choices))) {
      bp_choices <- character(0)
    }

    updateSelectizeInput(session, "bp1",
                         choices = bp_choices,
                         selected = NULL,
                         server = TRUE)
  })

  observeEvent(input$gen3, {
    req(input$gen3)

    bp_choices <- unique(datos()$breakpoint2[which(datos()$gene2 == input$gen3)])

    if (length(bp_choices) == 0 || all(is.na(bp_choices))) {
      bp_choices <- character(0)
    }

    updateSelectizeInput(session, "bp2",
                         choices = bp_choices,
                         selected = NULL,
                         server = TRUE)
  })


  fusion_results <- eventReactive(input$predict_btn, {
    req(input$gen5, input$gen3, input$bp1, input$bp2)

    # Convertir a entero
    bp1 <- as.integer(input$bp1)
    bp2 <- as.integer(input$bp2)

    SeqProt_Transcript <- Anotate_All_Transcripts(
      gene1 = input$gen5,
      gene2 = input$gen3,
      bp1   = bp1,
      bp2   = bp2,
      protein_df = protein_df
    )

    domains_df <- GetDomains_AllFusions(SeqProt_Transcript = SeqProt_Transcript)

    list(
      SeqProt_Transcript = SeqProt_Transcript,
      domains_df = domains_df
    )
  })

  # Mostrar las tablas
  output$tabla_fusion_transcript <- DT::renderDT({
    req(fusion_results())
    fusion_results()$SeqProt_Transcript
  }, options = list(scrollX = TRUE))

  output$tabla_fusion_domains <- DT::renderDT({
    req(fusion_results())
    df_list <- fusion_results()$domains_df

    # Verificar que sea lista y tenga al menos un elemento
    if (is.list(df_list) && length(df_list) >= 1) {
      df <- df_list[[1]]  # tomamos el primer data.frame
    } else {
      df <- data.frame(
        Info = paste0("Error: domains_df no es una lista con al menos un data.frame. Tipo recibido: ", class(df_list))
      )
    }

    df
  }, options = list(scrollX = TRUE))

  #--------------- QC -------------------------------

  observe({
    req(datos_qc())
    updateSelectInput(session, "selectCohorteQC",
                      choices = sort(unique(datos_qc()$Cohort)))
  })

  # Actualizar selectInput ID en función del Cohorte elegido
  observeEvent(input$selectCohorteQC, {
    req(datos_qc())
    ids <- datos_qc() %>% filter(Cohort == input$selectCohorteQC) %>% pull(ID) %>% unique()
    updateSelectInput(session, "selectIDQC", choices = sort(ids))
  })

  output$qc_secuenciacion_plot <- renderPlotly({
    req(input$selectCohorteQC, input$selectIDQC)

    # Filtrar por Cohorte y fuente QC Secuenciacion
    df <- datos_qc() %>%
      filter(Fuente == "QC Secuenciacion",
             Cohort == input$selectCohorteQC) %>%
      mutate(
        Valor = as.numeric(Valor),
        Metrica = factor(Metrica, levels = unique(Metrica)),
        color = ifelse(ID == input$selectIDQC, "red", "gray")  # marcar ID seleccionado
      ) %>%
      filter(!is.na(Valor))

    req(nrow(df) > 0)

    plot_ly(df, x = ~Metrica, y = ~Valor,
            type = "box",
            pointpos = 0,
            marker = list(color = ~color)) %>%
      layout(title = "QC Secuenciación",
             xaxis = list(title = "Métrica"),
             yaxis = list(title = "Valor"))
  })

  output$qc_alineamiento_plot <- renderPlotly({
    req(input$selectCohorteQC, input$selectIDQC)

    # Filtrar por Cohorte y Fuente
    df <- datos_qc() %>%
      filter(Fuente == "QC Alineamiento",
             Cohort == input$selectCohorteQC) %>%
      mutate(Valor = as.numeric(Valor))

    plots <- lapply(unique(df$Metrica), function(m) {
      df_m <- df %>% filter(Metrica == m)

      # Vector para colorear puntos: rojo si es el ID seleccionado, gris si no
      colors <- ifelse(df_m$ID == input$selectIDQC, "red", "gray")

      plot_ly(df_m, y = ~Valor, type = "box", name = m,
              boxpoints = df_m$ID[which(df_m$ID == input$selectIDQC)],
              #jitter = 0.3,
              pointpos = 0,
              marker = list(color = colors)) %>%
        layout(yaxis = list(title = m))
    })

    # Combinar subplots
    subplot(plots, nrows = 2, shareX = FALSE, shareY = FALSE, titleX = FALSE)
  })



}


# Ejecutar la aplicación
shinyApp(ui = ui, server = server)
