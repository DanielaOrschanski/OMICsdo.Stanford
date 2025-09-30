# FLENI Fusiones - Shiny App para Visualización y Anotación de Fusiones
# Cargar librerías necesarias
library(shiny)
library(shinydashboard)
library(DT)
library(readxl)
library(dplyr)
library(shinyFiles)
library(shinycssloaders)
library(plotly)

# Usar renderEchartsforR para que los graficos sean mas esteticos

# Ruta al archivo
ruta_archivo <- "/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/ReportFusions_TodasCohortes.xlsx"
datos <- read_excel(ruta_archivo)

# Ruta del archivo Excel
ruta_excel <- "/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/ReportFusions_TodasCohortes.xlsx"

# UI
ui <- dashboardPage(
  dashboardHeader(title = "FLENI - Análisis de Fusiones"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("by Patient", tabName = "fusiones", icon = icon("dna")),
      menuItem("by Fusion", tabName = "analisis", icon = icon("microscope"))
    )
  ),
  dashboardBody(
    tabItems(
      # --------------------- Dashboard ----------------------------------------------
      tabItem(tabName = "dashboard",
              fluidRow(
                box(width = 12,
                    title = "All Fusions",
                    status = "primary",
                    solidHeader = TRUE,
                    DTOutput("tabla")
                )
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
                    # Segunda fila de gráficos
                    fluidRow(
                      column(width = 6,
                             plotlyOutput("BoxplotCancer")
                      ),
                      column(width = 6,
                             plotlyOutput("BoxplotCancerConfidence")
                      )
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
              )
            ),
    # ---------------- by Patient ----------------
    tabItem(tabName = "fusiones",
            fluidRow(
              box(width = 12, title = "Patient Explorer", status = "primary", solidHeader = TRUE,

                  # Selectores dinámicos
                  fluidRow(
                    column(width = 6,
                           selectInput("selectCohorte", "Select Cohort:", choices = NULL)
                    ),
                    column(width = 6,
                           selectInput("selectID", "Select Patient ID:", choices = NULL)
                    )
                  ),

                  # Tabla de fusiones filtradas
                  DTOutput("tablaPaciente")
              )
            ),

            # Sub-ventanas cuando se selecciona una fusión
            fluidRow(
              box(width = 12, title = "Fusion details", status = "primary", solidHeader = TRUE,
                  tabBox(width = 12,
                         tabPanel("Fusion plot", plotOutput("fusionPlot")),

                         tabPanel("Domains",
                                  tabsetPanel(
                                    tabPanel("Gene 1 Domains", DT::dataTableOutput("domains_gene1")),
                                    tabPanel("Gene 2 Domains", DT::dataTableOutput("domains_gene2"))
                                  )
                         ),

                         tabPanel("Protein sequence", verbatimTextOutput("fusionSeq")),
                         #tabPanel("For docking", verbatimTextOutput("fusionDockingGene1")),
                         tabPanel("For docking",
                                  fluidPage(
                                    tabsetPanel(
                                      id = "fusionDockingTabs",
                                      type = "tabs",
                                      tabPanel("Gene1", verbatimTextOutput("fusionDockingGene1")),
                                      tabPanel("Gene2", verbatimTextOutput("fusionDockingGene2"))
                                    )
                                  )
                         ),

                         tabPanel("Drugs",
                                  tabsetPanel(
                                    tabPanel("Gene 1 Drugs", DT::dataTableOutput("drugs_gene1")),
                                    tabPanel("Gene 2 Drugs", DT::dataTableOutput("drugs_gene2"))
                                  )
                         ),
                  )
              )
            )


        )

    )
  )
)

#FUNCIONES A USAR EN EL SERVER --------------------------
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

      drugs_gene <- drugs_gene <- drugs_gene[1, which(colnames(drugs_gene) %in% c("gene_name", "drug_name"))]
      drugs_gene$drug_name <- "No drugs approved for this gene"

    } else {
      drugs_gene <- drugs_gene[, which(colnames(drugs_gene) %in% c("drug_name", "interaction_score", "interaction_type", "approved", "immunotherapy", "anti_neoplastic"))]
      if(any(duplicated(drugs_gene))) {
        drugs_gene <- drugs_gene[!duplicated(drugs_gene),]
      }
      drugs_gene <- drugs_gene[which(drugs_gene$approved == "TRUE"),]
      drugs_gene <- drugs_gene[, c(3,1,2,4,5,6)]
      rownames(drugs_gene) <- NULL
    }
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

  return(list(seq_prot_fusionada, sequence1, sequence2))
}

#transcript_id = "ENST00000646891"
#gene = "BRAF"
#sequence_gene = "MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEHIEALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTVTSSSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRDSLKKALMMRGLIPECCAVYRIQDGEKKPIGWDTDISWLTGEELHVEVLENVPLTTHNFVRKTFFTLAFCDFCRKLLFQGFRCQTCGYKFHQRCSTEVPLMCVNYDQLDLLFVSKFFEHHPIPQEEASLAETALTSGSSPSAPASDSIGPQILTSPSPSKSIPIPQPFRPADEDHRNQFGQRDRSSSAPNVHINTIEPVNIDDLIRDQGFRGDGGSTTGLSATPPASLPGSLTNVKALQKSPGPQRERKSSSSSEDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATEKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGAFPVH"
#seq_prot = "MAALSGGGGGGAEPGQALFNGDMEPEAGAGAGAAASSAADPAIPEEVWNIKQMIKLTQEHIEALLDKFGGEHNPPSIYLEAYEEYTSKLDALQQREQQLLESLGNGTDFSVSSSASMDTVTSSSSSSLSVLPSSLSVFQNPTDVARSNPKSPQKPIVRVFLPNKQRTVVPARCGVTVRDSLKKALMMRGLIPECCAVYRIQDGEKKPIGWDTDISWLTGEELHVEVLENVPLTTHNFVRKTFFTLAFCDFCRKLLFQGFRCQTCGYKFHQRCSTEVPLMCVNYDQLDLLFVSKFFEHHPIPQEEASLAETALTSGSSPSAPASDSIGPQILTSPSPSKSIPIPQPFRPADEDHRNQFGQRDRSSSAPNVHINTIEPVNIDDLIRDQGFRGDGGSTTGLSATPPASLPGSLTNVKALQKSPGPQRERKSSSSSEDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATEKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRSASEPSLNRAGFQTEDFSLYACASPKTPIQAGGYGAFPVH"

GetDomain <- function(transcript_id, gene, sequence_gene, seq_prot, path_dir = "/media/16TBDisk/Daniela/OMICsdo/R/MVP-Fusiones/fasta_tmp") {

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
    conserved_domain <- get_longest_common_substring(domain = domain, fusion = seq_prot)
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

# SERVER
server <- function(input, output, session) {

  # Leer el Excel al iniciar la app
  datos <- reactive({
    read_excel(ruta_excel)
  })

  # DASHBOARD ----------------------------------------------------------------

  # Mostrar tabla
  output$tabla <- renderDT({
    datatable(datos(), options = list(pageLength = 10, scrollX = TRUE))
  })

  # Gráfico de torta por cohorte (conteo de samples únicos)
  output$graficoCohorte <- renderPlotly({
    req(datos())

    resumen <- datos() %>%
      group_by(Cohorte) %>%
      summarise(SamplesUnicos = n_distinct(ID), .groups = "drop")

    plot_ly(
      data = resumen,
      labels = ~Cohorte,
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
      group_by(Cohorte) %>%
      summarise(Fusiones = n(), .groups = "drop")

    plot_ly(
      data = resumen_cohorte,
      x = ~Cohorte,
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
  output$BoxplotCancer <- renderPlotly({
    req(datos())

    resumen <- datos() %>%
      group_by(ID, Cancer) %>%
      summarise(FusionCount = n(), .groups = "drop")

    plot_ly(
      data = resumen,
      x = ~Cancer,
      y = ~FusionCount,
      type = "box"
      #boxpoints = "all"       # muestra puntos individuales
      #jitter = 0.5,            # dispersión de puntos
      #pointpos = -1.8          # posición de puntos
    ) %>%
      layout(title = "Fusions per ID by Cancer",
             xaxis = list(title = "Cancer"),
             yaxis = list(title = "Nº of Fusions"))
  })

  # Boxplot 2: Número de fusiones por ID, agrupado por Cancer y Confidence
  output$BoxplotCancerConfidence <- renderPlotly({
    req(datos())

    resumen <- datos() %>%
      group_by(ID, Cancer, confidence) %>%
      summarise(FusionCount = n(), .groups = "drop")

    plot_ly(
      data = resumen,
      x = ~Cancer,
      y = ~FusionCount,
      color = ~confidence,     # se colorea por Confidence
      type = "box"
      #boxpoints = "all"
      #jitter = 0.5,
      #pointpos = -1.8
    ) %>%
      layout(title = "Fusions per ID, Cancer and Confidence",
             xaxis = list(title = "Cancer"),
             yaxis = list(title = "Nº of Fusions"),
             boxmode = "confidence-group")   # agrupa por categoría de Confidence
  })

  #Barras con TK:
  output$TK_100_Cohorte <- renderPlotly({
    req(datos())

    df_counts <- datos() %>%
      mutate(Kinase = ifelse(grepl("Protein_kinase_domain\\(100%\\)", retained_protein_domains),
                             "Kinase domain", "Other")) %>%
      group_by(Cohorte, Kinase) %>%
      summarise(Fusions = n(), .groups = "drop")

    p <- ggplot(df_counts, aes(x = Cohorte, y = Fusions, fill = Kinase)) +
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
      group_by(Cohorte, Kinase) %>%
      summarise(Fusions = n(), .groups = "drop")

    p <- ggplot(df_counts, aes(x = Cohorte, y = Fusions, fill = Kinase)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = "Any % Retained Kinase Domain",
           x = "Cohort", y = "Fusion Count") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggplotly(p)
  })


  # BY PATIENT ---------------------------
  # Actualizar selectInput Cohorte
  observe({
    updateSelectInput(session, "selectCohorte",
                      choices = sort(unique(datos()$Cohorte)))
  })

  # Actualizar selectInput ID en función del Cohorte elegido
  observeEvent(input$selectCohorte, {
    ids <- datos() %>% filter(Cohorte == input$selectCohorte) %>% pull(ID) %>% unique()
    updateSelectInput(session, "selectID", choices = sort(ids))
  })

  # Filtrar tabla según Cohorte e ID
  datosPaciente <- reactive({
    req(input$selectCohorte, input$selectID)
    datos() %>% filter(Cohorte == input$selectCohorte, ID == input$selectID)
  })

  #output$tablaPaciente <- renderDT({
  #  datatable(datosPaciente(), selection = "single", options = list(scrollX = TRUE))
  #})
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
  output$fusionTabs <- renderUI({
    req(fusionesSeleccionadas())

    tabs <- lapply(1:nrow(fusionesSeleccionadas()), function(i) {
      fusion <- fusionesSeleccionadas()[i, ]
      fusion_id <- paste0("fusion_", i)  # id único

      tabPanel(
        title = paste("Fusion:", fusion$FusionName),
        tabsetPanel(
          tabPanel("Fusion plot", plotOutput(paste0("fusionPlot_", fusion_id))),
          tabPanel("Domains",
                   tabsetPanel(
                     tabPanel("Gene 1 Domains", DTOutput(paste0("domains_gene1_", fusion_id))),
                     tabPanel("Gene 2 Domains", DTOutput(paste0("domains_gene2_", fusion_id)))
                   )
          ),
          tabPanel("Protein sequence", verbatimTextOutput(paste0("fusionSeq_", fusion_id))),
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
      fusion_id <- paste0("fusion_", i)

      # Fusion plot
      local({
        f <- fusion
        id <- fusion_id
        output[[paste0("fusionPlot_", id)]] <- renderPlot({
          plot(1:10, main = paste("Fusion plot for", f$FusionName))
        })
      })

      # Domains Gene1
      local({
        f <- fusion
        id <- fusion_id
        output[[paste0("domains_gene1_", id)]] <- DT::renderDataTable({
          seq_total <- AnotateFusion(f$transcript_id1, f$transcript_id2, protein_df, f$peptide_sequence)
          seq_prot <- seq_total[[1]]
          sequence_gene1 <- seq_total[[2]]
          d <- GetDomain(f$transcript_id1, f$gene1, sequence_gene1, seq_prot)
          d[[1]]
        })
      })


    }
  })

  #############################
  # Subventana 1: gráfico de la fusión
  output$fusionPlot <- renderPlot({
    req(fusionSeleccionada())
    # Aquí usarías tu función propia, ej:
    # plotFusion(fusionSeleccionada())
    plot(1:10, main = paste("Fusion plot for", fusionSeleccionada()$FusionName))  # placeholder
  })

  # Subventana 2: dominios -----------------------------
  # Tabla de dominios gen1
  output$domains_gene1 <- DT::renderDataTable({
    req(fusionSeleccionada())
    withProgress(message = "Fusion domains are being calculated...", {
      seq_total <- AnotateFusion(
        transcript_id1 = fusionSeleccionada()$transcript_id1,
        transcript_id2 = fusionSeleccionada()$transcript_id2,
        protein_df     = protein_df,
        peptide_sequence = fusionSeleccionada()$peptide_sequence
      )

      seq_prot <- seq_total[[1]]
      sequence_gene1 <- seq_total[[2]]

      d <- GetDomain(
        transcript_id = fusionSeleccionada()$transcript_id1,
        gene = fusionSeleccionada()$gene1,
        sequence_gene = sequence_gene1,
        seq_prot = seq_prot
      )

      domains_df <- d[[1]]
      domains_df

    })
  }, options = list(
    scrollX = TRUE,   # barra horizontal
    scrollY = "400px", # barra vertical con altura fija
    paging   = FALSE   # sin paginación (opcional)
  ))


  # Tabla de dominios gen2
  output$domains_gene2 <- DT::renderDataTable({
    req(fusionSeleccionada())
    withProgress(message = "Fusion domains are being calculated...", {
      seq_total <- AnotateFusion(
        transcript_id1 = fusionSeleccionada()$transcript_id1,
        transcript_id2 = fusionSeleccionada()$transcript_id2,
        protein_df     = protein_df,
        peptide_sequence = fusionSeleccionada()$peptide_sequence
      )

      seq_prot <- seq_total[[1]]
      sequence_gene2 <- seq_total[[3]]

      d <- GetDomain(
        transcript_id = fusionSeleccionada()$transcript_id2,
        gene = fusionSeleccionada()$gene2,
        sequence_gene = sequence_gene2,
        seq_prot = seq_prot
      )
      domains_df <- d[[1]]
      domains_df

    })
  }, options = list(
    scrollX = TRUE,
    scrollY = "400px",
    paging   = FALSE
  ))


  # ------------------------------------------------

  # Subventana 3: secuencia proteína
  output$fusionSeq <- renderPrint({
    req(fusionSeleccionada())
    # Cargo RDS (ideal: hacerlo una sola vez en server y reusar)
    #BED_Anotador <- read_rds("/media/16TBDisk/Daniela/OMICsdo/data/BED_Anotador.RDS")
    protein_df   <- read_rds("/media/16TBDisk/Daniela/OMICsdo/data/DB_Anotador_Prot_DNA.RDS")

    seq_total <- AnotateFusion(
      transcript_id1= fusionSeleccionada()$transcript_id1,
      transcript_id2= fusionSeleccionada()$transcript_id2,
      protein_df    = protein_df,
      peptide_sequence = fusionSeleccionada()$peptide_sequence
    )

    seq_total
  })

  # Subventana 4: docking -----------------------------------
  # Gene1
  output$fusionDockingGene1 <- renderPrint({
    req(fusionSeleccionada())

    seq_total <- AnotateFusion(
      transcript_id1   = fusionSeleccionada()$transcript_id1,
      transcript_id2   = fusionSeleccionada()$transcript_id2,
      protein_df       = protein_df,
      peptide_sequence = fusionSeleccionada()$peptide_sequence
    )

    seq_prot       <- seq_total[[1]]
    sequence_gene1 <- seq_total[[2]]

    d <- GetDomain(
      transcript_id = fusionSeleccionada()$transcript_id1,
      gene          = fusionSeleccionada()$gene1,
      sequence_gene = sequence_gene1,
      seq_prot      = seq_prot
    )

    residuos <- d[[2]]
    residuos
  })


  # Gene2
  output$fusionDockingGene2 <- renderPrint({
    req(fusionSeleccionada())

    seq_total <- AnotateFusion(
      transcript_id1   = fusionSeleccionada()$transcript_id1,
      transcript_id2   = fusionSeleccionada()$transcript_id2,
      protein_df       = protein_df,
      peptide_sequence = fusionSeleccionada()$peptide_sequence
    )

    seq_prot       <- seq_total[[1]]
    sequence_gene2 <- seq_total[[3]]

    d <- GetDomain(
      transcript_id = fusionSeleccionada()$transcript_id2,
      gene          = fusionSeleccionada()$gene2,
      sequence_gene = sequence_gene2,
      seq_prot      = seq_prot
    )
    residuos <- d[[2]]
    residuos
  })


  # Subventana 5: Drogas ------------------------------------
  getDrugs(protein_df, gene)
  # Tabla de drogas gen1
  output$drugs_gene1 <- DT::renderDataTable({
    req(fusionSeleccionada())
    drugs_df1 <- getDrugs(protein_df = protein_df,
                          gene = fusionSeleccionada()$gene1)
  }, options = list(
    scrollY = TRUE,   # barra horizontal
    scrollX = "400px", # barra vertical con altura fija
    paging   = FALSE   # sin paginación (opcional)
  ))

  # Tabla de drogas gen2
  output$drugs_gene2 <- DT::renderDataTable({
    req(fusionSeleccionada())
    drugs_df2 <- getDrugs(protein_df = protein_df,
                          gene = fusionSeleccionada()$gene2)
  }, options = list(
    scrollY = TRUE,   # barra horizontal
    scrollX = "400px", # barra vertical con altura fija
    paging   = FALSE   # sin paginación (opcional)
  ))


}


# Ejecutar la aplicación
shinyApp(ui = ui, server = server)
