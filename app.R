## app.R ##
library(shiny)
library(shinydashboard)
library(tidyverse)
library(broom)
library(kableExtra)
library(readxl)
library(RCurl)
library(viridis)   

x <- getURL("https://raw.githubusercontent.com/demar01/rpmrna/master/S1_ProttobodyRNAseq.csv")
bp <- read.csv(text = x)
bp2<-bp %>%
select(-Gene.name) %>%
mutate_if(is.character,as.numeric)
bp2<-bp %>%
select(Gene.name ) %>%
cbind(bp2)
bp2<-bp2 %>% mutate_all(~replace(., is.na(.), 0))
bp3<- bp2 %>%
pivot_longer(-Gene.name, names_to="cell",values_to="log2") %>% 
    mutate_if(is.character,as.factor) 
    
bp3 %>% 
    filter(Gene.name=="AAAS") %>%
ggplot(aes(fct_inorder(cell),log2,fill= cell))+
    geom_col()+
    ylab("Log2 prot/body") + 
    xlab("") + 
    labs(title = "Log2 protrusion/cell body mRNA")+
    theme_minimal()+
    theme(legend.position = "none") 
    
proteins<- getURL("https://raw.githubusercontent.com/demar01/rpmrna/master/S2_ProttobodyTMT.csv")
proteins <- read.csv(text = proteins)
proteins<-proteins %>% mutate_all(~replace(., is.na(.), 0))
proteins3<- proteins %>%
    pivot_longer(-Protein, names_to="cell",values_to="log2") %>% 
    mutate_if(is.character,as.factor) %>% 
    filter(Protein!="") %>% arrange(Protein)
proteins3 %>% 
    filter(Protein=="LARP6") %>%
    ggplot(aes(fct_inorder(cell),log2,fill= cell))+
    geom_col()+
    ylab("Log2 prot/body") + 
    xlab("") + 
    labs(title = "Log2 (protrusion/cell body) Protein")+
    theme_minimal()+
    theme(legend.position = "none",
         text = element_text(size=20),
        axis.text.x = element_text(size=20)) 

ui <- dashboardPage(
    dashboardHeader(title = "Dermit et al. 2020"),
    
    dashboardSidebar(
        selectizeInput("v_gene", "mRNA",options = list(maxOptions = 10000),
                    choices = bp3 %>% select( Gene.name ) %>% distinct()),
        
        selectizeInput("v_protein", "Protein", options = list(maxOptions = 10000),
                       choices = 
                        proteins3 %>% select(Protein) %>% distinct())
    ),
    dashboardBody(
        fluidRow(box(plotOutput("prot_body")),box(plotOutput("prot_body_protein")))
    )
)

server <- function(input, output) { 
    
    output$prot_body <- renderPlot({
        bp3 %>% 
            filter(Gene.name==input$v_gene) %>%
            ggplot(aes(fct_inorder(cell),log2,fill= cell))+
            geom_col()+
            ylab("Log2 prot/body") + 
            xlab("") + 
            labs(title = "Log2 (protrusion/cell-body) mRNA")+
            theme_minimal()+
            theme(legend.position = "none") 
        
    })
    
    output$prot_body_protein <- renderPlot({
        proteins3 %>% 
            filter(Protein==input$v_protein) %>%
            ggplot(aes(fct_inorder(cell),log2,fill= cell))+
            geom_col()+
            ylab("Log2 prot/body") + 
            xlab("") + 
            labs(title = "Log2 (protrusion/cell-body) protein")+
            theme_minimal()+
            theme(legend.position = "none",
                  text = element_text(size=20),
        axis.text.x = element_text(size=20))
        
    })
    
    
    
}

shinyApp(ui, server)
