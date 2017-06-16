
#library(shinythemes)
library(shiny)
library(pipe.design)

shinyUI(fluidPage( # theme = shinytheme("cerulean") ,
  
  fluidRow(
    column(4 , titlePanel("PIPE Dose Escalation"))   ,
    column(8 , imageOutput("PIPEImage" , inline=TRUE , height="20%")) 
  ) ,
  
  tabsetPanel(
    tabPanel("Design",
             wellPanel( 
               fluidRow(
                 column(4 ,
                        numericInput("dimA" , value=2 , label="Dimension Drug A" , min=2 , max=8 , step=1) ,
                        numericInput("dimB" , value=2 , label="Dimension Drug B" , min=2 , max=8 , step=1) ,
                        numericInput("theta" , value=0.3 , label="Target DLT Rate" , min=0 , max=1 , step=0.05) ,
                        numericInput("cohortsize" , value=3 , label="Cohort Size" , min=1 , max=6 , step=1) 
                 ) ,
                 column(4 , 
                        selectInput("admis" , label="Admissible Doses" , choices=list("Closest Doses"="closest","Adjacent Doses"="adjacent") ,selected="closest") ,
                        selectInput("strategy" , label="Selecting amongst admissible doses" , choices=list("Minimum Sample Size"="ss","Weighted Randomisation (by sample size)"="ss-random"),selected="ss" ) ,
                        selectInput("constraint" , label="Dose Skipping Constraint" , choices=list("None"="none","Neighbouring Doses Only"="neighbouring","No Dose Skipping"="no.dose.skip","Neighbouring with no Diagonal Escalation"="neighbouring-nodiag","No Dose Skipping with no Diagonal Escalation"="no.dose.skip-nodiag"),selected="none") ,
                        helpText("`No Dose Skipping' allows next dose to be up to one dose level above any previously experimented drug A and drug B level"),
                        numericInput("epsilon" , value=1 , label="Safety Constraint (epsilon)" , min=0 , max=1 , step=0.05) ,
                        helpText("Doses whose weighted posterior probability of being greater than MTC is greater than epsilon will not be experimented on")          
                 ) ,
                 column(4 ,
                        selectInput("alternate" , label="Alternate" , choices=list("Yes"=T,"No"=F),selected=F) ,    
                        helpText("Should the trial always alternative above and below the MTC (subject to constraints)?"),
                        # selectInput("contourselect" , label="Contour Selection" , choices=c("mode","median")),
                        numericInput("seed" , value=1, label="Seed" , min=1, max=10000, step=1),
                        helpText("Starting seed to allow replication of the trial or simulation"),
                        numericInput("uppertox.constraint" , value=1 , label="Safety Constraint (upper toxicity)" , min=0 , max=1 , step=0.05) ,
                        helpText("No dose combination that lies above the most likely contour for this constraint may be dosed"),
                        numericInput("stop" , value=1 , label="Safety Constraint (stop)" , min=0 , max=1 , step=0.05) ,
                        helpText("Stop the trial if posterior probability of being greatrer than target DLT rate at lowest dose is higher than this number")
                 )
               ) 
             ), 
             wellPanel(
               fluidRow(
                 h4("Non-admissible doses")
               ),
               fluidRow(
                 column(4, helpText("Dose combinations that are non admissible throughout the trial (Drug A, Drug B)"))
               ),
               
               fluidRow(
                 column(1 , offset=1 ,
                        
                        conditionalPanel( condition = "input.dimA>=1 && input.dimB>=8" , checkboxInput("non.admissible_1_8" , value = F , label="(1,8)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=1 && input.dimB>=7" , checkboxInput("non.admissible_1_7" , value = F , label="(1,7)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=1 && input.dimB>=6" , checkboxInput("non.admissible_1_6" , value = F , label="(1,6)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=1 && input.dimB>=5" , checkboxInput("non.admissible_1_5" , value = F , label="(1,5)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=1 && input.dimB>=4" , checkboxInput("non.admissible_1_4" , value = F , label="(1,4)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=1 && input.dimB>=3" , checkboxInput("non.admissible_1_3" , value = F , label="(1,3)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=1 && input.dimB>=2" , checkboxInput("non.admissible_1_2" , value = F , label="(1,2)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=1 && input.dimB>=1" , checkboxInput("non.admissible_1_1" , value = F , label="(1,1)"  ) ) 
                 ) ,
                 column(1 , 
                        
                        conditionalPanel( condition = "input.dimA>=2 && input.dimB>=8" , checkboxInput("non.admissible_2_8" , value = F , label="(2,8)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=2 && input.dimB>=7" , checkboxInput("non.admissible_2_7" , value = F , label="(2,7)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=2 && input.dimB>=6" , checkboxInput("non.admissible_2_6" , value = F , label="(2,6)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=2 && input.dimB>=5" , checkboxInput("non.admissible_2_5" , value = F , label="(2,5)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=2 && input.dimB>=4" , checkboxInput("non.admissible_2_4" , value = F , label="(2,4)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=2 && input.dimB>=3" , checkboxInput("non.admissible_2_3" , value = F , label="(2,3)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=2 && input.dimB>=2" , checkboxInput("non.admissible_2_2" , value = F , label="(2,2)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=2 && input.dimB>=1" , checkboxInput("non.admissible_2_1" , value = F , label="(2,1)"  ) ) 
                 )  ,  
                 column(1 , 
                        
                        conditionalPanel( condition = "input.dimA>=3 && input.dimB>=8" , checkboxInput("non.admissible_3_8" , value = F , label="(3,8)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=3 && input.dimB>=7" , checkboxInput("non.admissible_3_7" , value = F , label="(3,7)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=3 && input.dimB>=6" , checkboxInput("non.admissible_3_6" , value = F , label="(3,6)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=3 && input.dimB>=5" , checkboxInput("non.admissible_3_5" , value = F , label="(3,5)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=3 && input.dimB>=4" , checkboxInput("non.admissible_3_4" , value = F , label="(3,4)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=3 && input.dimB>=3" , checkboxInput("non.admissible_3_3" , value = F , label="(3,3)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=3 && input.dimB>=2" , checkboxInput("non.admissible_3_2" , value = F , label="(3,2)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=3 && input.dimB>=1" , checkboxInput("non.admissible_3_1" , value = F , label="(3,1)"  ) ) 
                 ) ,
                 column(1 , 
                        
                        conditionalPanel( condition = "input.dimA>=4 && input.dimB>=8" , checkboxInput("non.admissible_4_8" , value = F , label="(4,8)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=4 && input.dimB>=7" , checkboxInput("non.admissible_4_7" , value = F , label="(4,7)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=4 && input.dimB>=6" , checkboxInput("non.admissible_4_6" , value = F , label="(4,6)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=4 && input.dimB>=5" , checkboxInput("non.admissible_4_5" , value = F , label="(4,5)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=4 && input.dimB>=4" , checkboxInput("non.admissible_4_4" , value = F , label="(4,4)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=4 && input.dimB>=3" , checkboxInput("non.admissible_4_3" , value = F , label="(4,3)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=4 && input.dimB>=2" , checkboxInput("non.admissible_4_2" , value = F , label="(4,2)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=4 && input.dimB>=1" , checkboxInput("non.admissible_4_1" , value = F , label="(4,1)"  ) ) 
                 )  ,  
                 column(1 , 
                        
                        conditionalPanel( condition = "input.dimA>=5 && input.dimB>=8" , checkboxInput("non.admissible_5_8" , value = F , label="(5,8)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=5 && input.dimB>=7" , checkboxInput("non.admissible_5_7" , value = F , label="(5,7)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=5 && input.dimB>=6" , checkboxInput("non.admissible_5_6" , value = F , label="(5,6)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=5 && input.dimB>=5" , checkboxInput("non.admissible_5_5" , value = F , label="(5,5)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=5 && input.dimB>=4" , checkboxInput("non.admissible_5_4" , value = F , label="(5,4)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=5 && input.dimB>=3" , checkboxInput("non.admissible_5_3" , value = F , label="(5,3)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=5 && input.dimB>=2" , checkboxInput("non.admissible_5_2" , value = F , label="(5,2)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=5 && input.dimB>=1" , checkboxInput("non.admissible_5_1" , value = F , label="(5,1)"  ) ) 
                 ) ,
                 column(1 , 
                        
                        conditionalPanel( condition = "input.dimA>=6 && input.dimB>=8" , checkboxInput("non.admissible_6_8" , value = F , label="(6,8)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=6 && input.dimB>=7" , checkboxInput("non.admissible_6_7" , value = F , label="(6,7)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=6 && input.dimB>=6" , checkboxInput("non.admissible_6_6" , value = F , label="(6,6)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=6 && input.dimB>=5" , checkboxInput("non.admissible_6_5" , value = F , label="(6,5)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=6 && input.dimB>=4" , checkboxInput("non.admissible_6_4" , value = F , label="(6,4)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=6 && input.dimB>=3" , checkboxInput("non.admissible_6_3" , value = F , label="(6,3)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=6 && input.dimB>=2" , checkboxInput("non.admissible_6_2" , value = F , label="(6,2)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=6 && input.dimB>=1" , checkboxInput("non.admissible_6_1" , value = F , label="(6,1)"  ) ) 
                 )  ,  
                 column(1 , 
                        conditionalPanel( condition = "input.dimA>=7 && input.dimB>=8" , checkboxInput("non.admissible_7_8" , value = F , label="(7,8)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=7 && input.dimB>=7" , checkboxInput("non.admissible_7_7" , value = F , label="(7,7)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=7 && input.dimB>=6" , checkboxInput("non.admissible_7_6" , value = F , label="(7,6)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=7 && input.dimB>=5" , checkboxInput("non.admissible_7_5" , value = F , label="(7,5)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=7 && input.dimB>=4" , checkboxInput("non.admissible_7_4" , value = F , label="(7,4)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=7 && input.dimB>=3" , checkboxInput("non.admissible_7_3" , value = F , label="(7,3)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=7 && input.dimB>=2" , checkboxInput("non.admissible_7_2" , value = F , label="(7,2)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=7 && input.dimB>=1" , checkboxInput("non.admissible_7_1" , value = F , label="(7,1)"  ) ) 
                 ) ,
                 column(1 , 
                        conditionalPanel( condition = "input.dimA>=8 && input.dimB>=8" , checkboxInput("non.admissible_8_8" , value = F , label="(8,8)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=8 && input.dimB>=7" , checkboxInput("non.admissible_8_7" , value = F , label="(8,7)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=8 && input.dimB>=6" , checkboxInput("non.admissible_8_6" , value = F , label="(8,6)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=8 && input.dimB>=5" , checkboxInput("non.admissible_8_5" , value = F , label="(8,5)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=8 && input.dimB>=4" , checkboxInput("non.admissible_8_4" , value = F , label="(8,4)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=8 && input.dimB>=3" , checkboxInput("non.admissible_8_3" , value = F , label="(8,3)"  ) ) ,
                        conditionalPanel( condition = "input.dimA>=8 && input.dimB>=2" , checkboxInput("non.admissible_8_2" , value = F , label="(8,2)"  ) ) ,                              
                        conditionalPanel( condition = "input.dimA>=8 && input.dimB>=1" , checkboxInput("non.admissible_8_1" , value = F , label="(8,1)"  ) ) 
                 ) 
               ))             
             
             
    ),
    
    tabPanel("Priors / Historical data",
             
             checkboxInput("pmedian", "Specify prior medians and sample sizes"),
             checkboxInput("p_a_and_b","Specify prior DLTs / no DLTs"),
             conditionalPanel(
               condition = "input.pmedian == true",
               wellPanel(
                 fluidRow(
                   h4("Prior medians")
                 ),
                 fluidRow(
                   column(3 , offset=3 , actionButton("interpolate","Interpolate Prior Medians",icon=icon("arrows-alt"))) ,
                   column(1 , offset=0 , actionButton("flat.prior","Flat Priors"))
                 ),
                 
                 fluidRow(
                   column(4, helpText("Prior Medians (Drug A, Drug B)"))
                 ),
                 
                 fluidRow(
                   column(1 , offset=1 ,
                          
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=8" , numericInput("prior_1_8" , value=NA , label="(1,8)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=7" , numericInput("prior_1_7" , value=NA , label="(1,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=6" , numericInput("prior_1_6" , value=NA , label="(1,6)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=5" , numericInput("prior_1_5" , value=NA , label="(1,5)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=4" , numericInput("prior_1_4" , value=NA , label="(1,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=3" , numericInput("prior_1_3" , value=NA , label="(1,3)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=2" , numericInput("prior_1_2" , value=NA , label="(1,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=1" , numericInput("prior_1_1" , value=NA , label="(1,1)" , min=0 , max=1 , step=0.05 ) ) 
                   ) ,
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=8" , numericInput("prior_2_8" , value=NA , label="(2,8)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=7" , numericInput("prior_2_7" , value=NA , label="(2,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=6" , numericInput("prior_2_6" , value=NA , label="(2,6)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=5" , numericInput("prior_2_5" , value=NA , label="(2,5)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=4" , numericInput("prior_2_4" , value=NA , label="(2,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=3" , numericInput("prior_2_3" , value=NA , label="(2,3)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=2" , numericInput("prior_2_2" , value=NA , label="(2,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=1" , numericInput("prior_2_1" , value=NA , label="(2,1)" , min=0 , max=1 , step=0.05 ) ) 
                   )  ,  
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=8" , numericInput("prior_3_8" , value=NA , label="(3,8)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=7" , numericInput("prior_3_7" , value=NA , label="(3,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=6" , numericInput("prior_3_6" , value=NA , label="(3,6)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=5" , numericInput("prior_3_5" , value=NA , label="(3,5)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=4" , numericInput("prior_3_4" , value=NA , label="(3,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=3" , numericInput("prior_3_3" , value=NA , label="(3,3)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=2" , numericInput("prior_3_2" , value=NA , label="(3,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=1" , numericInput("prior_3_1" , value=NA , label="(3,1)" , min=0 , max=1 , step=0.05 ) ) 
                   ) ,
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=8" , numericInput("prior_4_8" , value=NA , label="(4,8)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=7" , numericInput("prior_4_7" , value=NA , label="(4,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=6" , numericInput("prior_4_6" , value=NA , label="(4,6)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=5" , numericInput("prior_4_5" , value=NA , label="(4,5)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=4" , numericInput("prior_4_4" , value=NA , label="(4,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=3" , numericInput("prior_4_3" , value=NA , label="(4,3)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=2" , numericInput("prior_4_2" , value=NA , label="(4,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=1" , numericInput("prior_4_1" , value=NA , label="(4,1)" , min=0 , max=1 , step=0.05 ) ) 
                   )  ,  
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=8" , numericInput("prior_5_8" , value=NA , label="(5,8)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=7" , numericInput("prior_5_7" , value=NA , label="(5,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=6" , numericInput("prior_5_6" , value=NA , label="(5,6)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=5" , numericInput("prior_5_5" , value=NA , label="(5,5)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=4" , numericInput("prior_5_4" , value=NA , label="(5,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=3" , numericInput("prior_5_3" , value=NA , label="(5,3)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=2" , numericInput("prior_5_2" , value=NA , label="(5,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=1" , numericInput("prior_5_1" , value=NA , label="(5,1)" , min=0 , max=1 , step=0.05 ) ) 
                   ) ,
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=8" , numericInput("prior_6_8" , value=NA , label="(6,8)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=7" , numericInput("prior_6_7" , value=NA , label="(6,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=6" , numericInput("prior_6_6" , value=NA , label="(6,6)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=5" , numericInput("prior_6_5" , value=NA , label="(6,5)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=4" , numericInput("prior_6_4" , value=NA , label="(6,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=3" , numericInput("prior_6_3" , value=NA , label="(6,3)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=2" , numericInput("prior_6_2" , value=NA , label="(6,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=1" , numericInput("prior_6_1" , value=NA , label="(6,1)" , min=0 , max=1 , step=0.05 ) ) 
                   )  ,  
                   column(1 , 
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=8" , numericInput("prior_7_8" , value=NA , label="(7,8)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=7" , numericInput("prior_7_7" , value=NA , label="(7,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=6" , numericInput("prior_7_6" , value=NA , label="(7,6)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=5" , numericInput("prior_7_5" , value=NA , label="(7,5)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=4" , numericInput("prior_7_4" , value=NA , label="(7,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=3" , numericInput("prior_7_3" , value=NA , label="(7,3)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=2" , numericInput("prior_7_2" , value=NA , label="(7,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=1" , numericInput("prior_7_1" , value=NA , label="(7,1)" , min=0 , max=1 , step=0.05 ) ) 
                   ) ,
                   column(1 , 
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=8" , numericInput("prior_8_8" , value=NA , label="(8,8)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=7" , numericInput("prior_8_7" , value=NA , label="(8,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=6" , numericInput("prior_8_6" , value=NA , label="(8,6)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=5" , numericInput("prior_8_5" , value=NA , label="(8,5)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=4" , numericInput("prior_8_4" , value=NA , label="(8,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=3" , numericInput("prior_8_3" , value=NA , label="(8,3)" , min=0 , max=1 , step=0.05 ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=2" , numericInput("prior_8_2" , value=NA , label="(8,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=1" , numericInput("prior_8_1" , value=NA , label="(8,1)" , min=0 , max=1 , step=0.05 ) ) 
                   ) 
                 )),
               
               
               ## Prior Sample Sizes
               wellPanel(
                 fluidRow(
                   h4("Prior sample sizes")
                 ),
                 
                 fluidRow(
                   #column(3 , offset=3 , actionButton("totalss","Use same sample size for all dose combinations")) 
                   checkboxInput("totalss", "Use same sample size for all dose combinations")
                 ),
                 
                 conditionalPanel(
                   condition = "input.totalss==true",
                   numericInput("priorss" , value=1 , label="Total Prior Sample Size (over all dose combinations)" , min=0),
                   helpText("For a weak prior set the total prior sample size to 1, whilst for a strong prior set the total prior sample size equal to the number of dose combinations")
                 ),
                 
                 fluidRow(
                   column(4, helpText("Prior Sample Sizes (Drug A, Drug B)"))
                 ),
                 
                 fluidRow(
                   column(1 , offset=1 ,
                          
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=8" , numericInput("priorss_1_8" , value=NA , label="(1,8)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=7" , numericInput("priorss_1_7" , value=NA , label="(1,7)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=6" , numericInput("priorss_1_6" , value=NA , label="(1,6)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=5" , numericInput("priorss_1_5" , value=NA , label="(1,5)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=4" , numericInput("priorss_1_4" , value=NA , label="(1,4)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=3" , numericInput("priorss_1_3" , value=NA , label="(1,3)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=2" , numericInput("priorss_1_2" , value=NA , label="(1,2)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=1" , numericInput("priorss_1_1" , value=NA , label="(1,1)" , min=1e-8   ) ) 
                   ) ,
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=8" , numericInput("priorss_2_8" , value=NA , label="(2,8)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=7" , numericInput("priorss_2_7" , value=NA , label="(2,7)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=6" , numericInput("priorss_2_6" , value=NA , label="(2,6)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=5" , numericInput("priorss_2_5" , value=NA , label="(2,5)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=4" , numericInput("priorss_2_4" , value=NA , label="(2,4)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=3" , numericInput("priorss_2_3" , value=NA , label="(2,3)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=2" , numericInput("priorss_2_2" , value=NA , label="(2,2)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=1" , numericInput("priorss_2_1" , value=NA , label="(2,1)" , min=1e-8   ) ) 
                   )  ,  
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=8" , numericInput("priorss_3_8" , value=NA , label="(3,8)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=7" , numericInput("priorss_3_7" , value=NA , label="(3,7)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=6" , numericInput("priorss_3_6" , value=NA , label="(3,6)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=5" , numericInput("priorss_3_5" , value=NA , label="(3,5)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=4" , numericInput("priorss_3_4" , value=NA , label="(3,4)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=3" , numericInput("priorss_3_3" , value=NA , label="(3,3)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=2" , numericInput("priorss_3_2" , value=NA , label="(3,2)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=1" , numericInput("priorss_3_1" , value=NA , label="(3,1)" , min=1e-8   ) ) 
                   ) ,
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=8" , numericInput("priorss_4_8" , value=NA , label="(4,8)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=7" , numericInput("priorss_4_7" , value=NA , label="(4,7)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=6" , numericInput("priorss_4_6" , value=NA , label="(4,6)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=5" , numericInput("priorss_4_5" , value=NA , label="(4,5)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=4" , numericInput("priorss_4_4" , value=NA , label="(4,4)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=3" , numericInput("priorss_4_3" , value=NA , label="(4,3)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=2" , numericInput("priorss_4_2" , value=NA , label="(4,2)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=1" , numericInput("priorss_4_1" , value=NA , label="(4,1)" , min=1e-8   ) ) 
                   )  ,  
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=8" , numericInput("priorss_5_8" , value=NA , label="(5,8)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=7" , numericInput("priorss_5_7" , value=NA , label="(5,7)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=6" , numericInput("priorss_5_6" , value=NA , label="(5,6)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=5" , numericInput("priorss_5_5" , value=NA , label="(5,5)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=4" , numericInput("priorss_5_4" , value=NA , label="(5,4)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=3" , numericInput("priorss_5_3" , value=NA , label="(5,3)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=2" , numericInput("priorss_5_2" , value=NA , label="(5,2)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=1" , numericInput("priorss_5_1" , value=NA , label="(5,1)" , min=1e-8   ) ) 
                   ) ,
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=8" , numericInput("priorss_6_8" , value=NA , label="(6,8)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=7" , numericInput("priorss_6_7" , value=NA , label="(6,7)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=6" , numericInput("priorss_6_6" , value=NA , label="(6,6)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=5" , numericInput("priorss_6_5" , value=NA , label="(6,5)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=4" , numericInput("priorss_6_4" , value=NA , label="(6,4)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=3" , numericInput("priorss_6_3" , value=NA , label="(6,3)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=2" , numericInput("priorss_6_2" , value=NA , label="(6,2)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=1" , numericInput("priorss_6_1" , value=NA , label="(6,1)" , min=1e-8   ) ) 
                   )  ,  
                   column(1 , 
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=8" , numericInput("priorss_7_8" , value=NA , label="(7,8)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=7" , numericInput("priorss_7_7" , value=NA , label="(7,7)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=6" , numericInput("priorss_7_6" , value=NA , label="(7,6)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=5" , numericInput("priorss_7_5" , value=NA , label="(7,5)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=4" , numericInput("priorss_7_4" , value=NA , label="(7,4)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=3" , numericInput("priorss_7_3" , value=NA , label="(7,3)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=2" , numericInput("priorss_7_2" , value=NA , label="(7,2)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=1" , numericInput("priorss_7_1" , value=NA , label="(7,1)" , min=1e-8   ) ) 
                   ) ,
                   column(1 , 
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=8" , numericInput("priorss_8_8" , value=NA , label="(8,8)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=7" , numericInput("priorss_8_7" , value=NA , label="(8,7)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=6" , numericInput("priorss_8_6" , value=NA , label="(8,6)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=5" , numericInput("priorss_8_5" , value=NA , label="(8,5)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=4" , numericInput("priorss_8_4" , value=NA , label="(8,4)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=3" , numericInput("priorss_8_3" , value=NA , label="(8,3)" , min=1e-8   ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=2" , numericInput("priorss_8_2" , value=NA , label="(8,2)" , min=1e-8   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=1" , numericInput("priorss_8_1" , value=NA , label="(8,1)" , min=1e-8   ) ) 
                   ) 
                 ))),
             conditionalPanel(
               condition = "input.p_a_and_b == true",
               wellPanel(
                 fluidRow(
                   h4("Prior number of DLTs")
                 ),

                 fluidRow(
                   column(4, helpText("Prior number of DLTs (Drug A, Drug B)"))
                 ),
                 
                 fluidRow(
                   column(1 , offset=1 ,
                          
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=8" , numericInput("priora_1_8" , value=NA , label="(1,8)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=7" , numericInput("priora_1_7" , value=NA , label="(1,7)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=6" , numericInput("priora_1_6" , value=NA , label="(1,6)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=5" , numericInput("priora_1_5" , value=NA , label="(1,5)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=4" , numericInput("priora_1_4" , value=NA , label="(1,4)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=3" , numericInput("priora_1_3" , value=NA , label="(1,3)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=2" , numericInput("priora_1_2" , value=NA , label="(1,2)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=1" , numericInput("priora_1_1" , value=NA , label="(1,1)" , min=0 , max=999 , step=1 ) ) 
                   ) ,
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=8" , numericInput("priora_2_8" , value=NA , label="(2,8)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=7" , numericInput("priora_2_7" , value=NA , label="(2,7)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=6" , numericInput("priora_2_6" , value=NA , label="(2,6)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=5" , numericInput("priora_2_5" , value=NA , label="(2,5)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=4" , numericInput("priora_2_4" , value=NA , label="(2,4)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=3" , numericInput("priora_2_3" , value=NA , label="(2,3)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=2" , numericInput("priora_2_2" , value=NA , label="(2,2)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=1" , numericInput("priora_2_1" , value=NA , label="(2,1)" , min=0 , max=999 , step=1 ) ) 
                   )  ,  
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=8" , numericInput("priora_3_8" , value=NA , label="(3,8)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=7" , numericInput("priora_3_7" , value=NA , label="(3,7)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=6" , numericInput("priora_3_6" , value=NA , label="(3,6)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=5" , numericInput("priora_3_5" , value=NA , label="(3,5)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=4" , numericInput("priora_3_4" , value=NA , label="(3,4)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=3" , numericInput("priora_3_3" , value=NA , label="(3,3)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=2" , numericInput("priora_3_2" , value=NA , label="(3,2)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=1" , numericInput("priora_3_1" , value=NA , label="(3,1)" , min=0 , max=999 , step=1 ) ) 
                   ) ,
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=8" , numericInput("priora_4_8" , value=NA , label="(4,8)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=7" , numericInput("priora_4_7" , value=NA , label="(4,7)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=6" , numericInput("priora_4_6" , value=NA , label="(4,6)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=5" , numericInput("priora_4_5" , value=NA , label="(4,5)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=4" , numericInput("priora_4_4" , value=NA , label="(4,4)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=3" , numericInput("priora_4_3" , value=NA , label="(4,3)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=2" , numericInput("priora_4_2" , value=NA , label="(4,2)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=1" , numericInput("priora_4_1" , value=NA , label="(4,1)" , min=0 , max=999 , step=1 ) ) 
                   )  ,  
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=8" , numericInput("priora_5_8" , value=NA , label="(5,8)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=7" , numericInput("priora_5_7" , value=NA , label="(5,7)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=6" , numericInput("priora_5_6" , value=NA , label="(5,6)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=5" , numericInput("priora_5_5" , value=NA , label="(5,5)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=4" , numericInput("priora_5_4" , value=NA , label="(5,4)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=3" , numericInput("priora_5_3" , value=NA , label="(5,3)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=2" , numericInput("priora_5_2" , value=NA , label="(5,2)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=1" , numericInput("priora_5_1" , value=NA , label="(5,1)" , min=0 , max=999 , step=1 ) ) 
                   ) ,
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=8" , numericInput("priora_6_8" , value=NA , label="(6,8)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=7" , numericInput("priora_6_7" , value=NA , label="(6,7)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=6" , numericInput("priora_6_6" , value=NA , label="(6,6)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=5" , numericInput("priora_6_5" , value=NA , label="(6,5)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=4" , numericInput("priora_6_4" , value=NA , label="(6,4)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=3" , numericInput("priora_6_3" , value=NA , label="(6,3)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=2" , numericInput("priora_6_2" , value=NA , label="(6,2)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=1" , numericInput("priora_6_1" , value=NA , label="(6,1)" , min=0 , max=999 , step=1 ) ) 
                   )  ,  
                   column(1 , 
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=8" , numericInput("priora_7_8" , value=NA , label="(7,8)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=7" , numericInput("priora_7_7" , value=NA , label="(7,7)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=6" , numericInput("priora_7_6" , value=NA , label="(7,6)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=5" , numericInput("priora_7_5" , value=NA , label="(7,5)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=4" , numericInput("priora_7_4" , value=NA , label="(7,4)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=3" , numericInput("priora_7_3" , value=NA , label="(7,3)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=2" , numericInput("priora_7_2" , value=NA , label="(7,2)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=1" , numericInput("priora_7_1" , value=NA , label="(7,1)" , min=0 , max=999 , step=1 ) ) 
                   ) ,
                   column(1 , 
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=8" , numericInput("priora_8_8" , value=NA , label="(8,8)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=7" , numericInput("priora_8_7" , value=NA , label="(8,7)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=6" , numericInput("priora_8_6" , value=NA , label="(8,6)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=5" , numericInput("priora_8_5" , value=NA , label="(8,5)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=4" , numericInput("priora_8_4" , value=NA , label="(8,4)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=3" , numericInput("priora_8_3" , value=NA , label="(8,3)" , min=0 , max=999 , step=1 ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=2" , numericInput("priora_8_2" , value=NA , label="(8,2)" , min=0 , max=999 , step=1 ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=1" , numericInput("priora_8_1" , value=NA , label="(8,1)" , min=0 , max=999 , step=1 ) ) 
                   ) 
                 )),
               
               
               ## Prior No DLTs
               wellPanel(
                 fluidRow(
                   h4("Prior number of  No DLTs")
                 ),
                 
                 
                 fluidRow(
                   column(4, helpText("Prior number of No DLTs (Drug A, Drug B)"))
                 ),
                 
                 fluidRow(
                   column(1 , offset=1 ,
                          
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=8" , numericInput("priorb_1_8" , value=NA , label="(1,8)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=7" , numericInput("priorb_1_7" , value=NA , label="(1,7)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=6" , numericInput("priorb_1_6" , value=NA , label="(1,6)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=5" , numericInput("priorb_1_5" , value=NA , label="(1,5)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=4" , numericInput("priorb_1_4" , value=NA , label="(1,4)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=3" , numericInput("priorb_1_3" , value=NA , label="(1,3)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=2" , numericInput("priorb_1_2" , value=NA , label="(1,2)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=1 && input.dimB>=1" , numericInput("priorb_1_1" , value=NA , label="(1,1)" , min=0 , max=999 , step=1   ) ) 
                   ) ,
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=8" , numericInput("priorb_2_8" , value=NA , label="(2,8)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=7" , numericInput("priorb_2_7" , value=NA , label="(2,7)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=6" , numericInput("priorb_2_6" , value=NA , label="(2,6)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=5" , numericInput("priorb_2_5" , value=NA , label="(2,5)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=4" , numericInput("priorb_2_4" , value=NA , label="(2,4)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=3" , numericInput("priorb_2_3" , value=NA , label="(2,3)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=2" , numericInput("priorb_2_2" , value=NA , label="(2,2)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=2 && input.dimB>=1" , numericInput("priorb_2_1" , value=NA , label="(2,1)" , min=0 , max=999 , step=1   ) ) 
                   )  ,  
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=8" , numericInput("priorb_3_8" , value=NA , label="(3,8)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=7" , numericInput("priorb_3_7" , value=NA , label="(3,7)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=6" , numericInput("priorb_3_6" , value=NA , label="(3,6)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=5" , numericInput("priorb_3_5" , value=NA , label="(3,5)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=4" , numericInput("priorb_3_4" , value=NA , label="(3,4)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=3" , numericInput("priorb_3_3" , value=NA , label="(3,3)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=2" , numericInput("priorb_3_2" , value=NA , label="(3,2)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=3 && input.dimB>=1" , numericInput("priorb_3_1" , value=NA , label="(3,1)" , min=0 , max=999 , step=1   ) ) 
                   ) ,
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=8" , numericInput("priorb_4_8" , value=NA , label="(4,8)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=7" , numericInput("priorb_4_7" , value=NA , label="(4,7)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=6" , numericInput("priorb_4_6" , value=NA , label="(4,6)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=5" , numericInput("priorb_4_5" , value=NA , label="(4,5)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=4" , numericInput("priorb_4_4" , value=NA , label="(4,4)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=3" , numericInput("priorb_4_3" , value=NA , label="(4,3)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=2" , numericInput("priorb_4_2" , value=NA , label="(4,2)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=4 && input.dimB>=1" , numericInput("priorb_4_1" , value=NA , label="(4,1)" , min=0 , max=999 , step=1   ) ) 
                   )  ,  
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=8" , numericInput("priorb_5_8" , value=NA , label="(5,8)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=7" , numericInput("priorb_5_7" , value=NA , label="(5,7)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=6" , numericInput("priorb_5_6" , value=NA , label="(5,6)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=5" , numericInput("priorb_5_5" , value=NA , label="(5,5)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=4" , numericInput("priorb_5_4" , value=NA , label="(5,4)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=3" , numericInput("priorb_5_3" , value=NA , label="(5,3)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=2" , numericInput("priorb_5_2" , value=NA , label="(5,2)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=5 && input.dimB>=1" , numericInput("priorb_5_1" , value=NA , label="(5,1)" , min=0 , max=999 , step=1   ) ) 
                   ) ,
                   column(1 , 
                          
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=8" , numericInput("priorb_6_8" , value=NA , label="(6,8)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=7" , numericInput("priorb_6_7" , value=NA , label="(6,7)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=6" , numericInput("priorb_6_6" , value=NA , label="(6,6)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=5" , numericInput("priorb_6_5" , value=NA , label="(6,5)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=4" , numericInput("priorb_6_4" , value=NA , label="(6,4)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=3" , numericInput("priorb_6_3" , value=NA , label="(6,3)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=2" , numericInput("priorb_6_2" , value=NA , label="(6,2)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=6 && input.dimB>=1" , numericInput("priorb_6_1" , value=NA , label="(6,1)" , min=0 , max=999 , step=1   ) ) 
                   )  ,  
                   column(1 , 
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=8" , numericInput("priorb_7_8" , value=NA , label="(7,8)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=7" , numericInput("priorb_7_7" , value=NA , label="(7,7)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=6" , numericInput("priorb_7_6" , value=NA , label="(7,6)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=5" , numericInput("priorb_7_5" , value=NA , label="(7,5)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=4" , numericInput("priorb_7_4" , value=NA , label="(7,4)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=3" , numericInput("priorb_7_3" , value=NA , label="(7,3)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=2" , numericInput("priorb_7_2" , value=NA , label="(7,2)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=7 && input.dimB>=1" , numericInput("priorb_7_1" , value=NA , label="(7,1)" , min=0 , max=999 , step=1   ) ) 
                   ) ,
                   column(1 , 
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=8" , numericInput("priorb_8_8" , value=NA , label="(8,8)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=7" , numericInput("priorb_8_7" , value=NA , label="(8,7)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=6" , numericInput("priorb_8_6" , value=NA , label="(8,6)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=5" , numericInput("priorb_8_5" , value=NA , label="(8,5)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=4" , numericInput("priorb_8_4" , value=NA , label="(8,4)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=3" , numericInput("priorb_8_3" , value=NA , label="(8,3)" , min=0 , max=999 , step=1   ) ) ,
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=2" , numericInput("priorb_8_2" , value=NA , label="(8,2)" , min=0 , max=999 , step=1   ) ) ,                              
                          conditionalPanel( condition = "input.dimA>=8 && input.dimB>=1" , numericInput("priorb_8_1" , value=NA , label="(8,1)" , min=0 , max=999 , step=1   ) ) 
                   ) 
                 )))
             
    ),
    tabPanel("Run Trial",
             fluidRow(
               column(1 , offset=0 , actionButton("reset","Reset Study",icon=icon("repeat")))    
             ) ,
             
             fluidRow(
               column(1 , 
                      conditionalPanel( condition = "output.n>=1" ,uiOutput("doseA_1")) ,
                      conditionalPanel( condition = "output.n>=2" ,uiOutput("doseA_2")) ,
                      conditionalPanel( condition = "output.n>=3" ,uiOutput("doseA_3")) ,
                      conditionalPanel( condition = "output.n>=4" ,uiOutput("doseA_4")) ,
                      conditionalPanel( condition = "output.n>=5" ,uiOutput("doseA_5")) ,
                      conditionalPanel( condition = "output.n>=6" ,uiOutput("doseA_6")) ,
                      conditionalPanel( condition = "output.n>=7" ,uiOutput("doseA_7")) ,
                      conditionalPanel( condition = "output.n>=8" ,uiOutput("doseA_8")) ,
                      conditionalPanel( condition = "output.n>=9" ,uiOutput("doseA_9")) ,
                      conditionalPanel( condition = "output.n>=10" ,uiOutput("doseA_10")) ,
                      conditionalPanel( condition = "output.n>=11" ,uiOutput("doseA_11")) ,
                      conditionalPanel( condition = "output.n>=12" ,uiOutput("doseA_12")) ,
                      conditionalPanel( condition = "output.n>=13" ,uiOutput("doseA_13")) ,
                      conditionalPanel( condition = "output.n>=14" ,uiOutput("doseA_14")) ,
                      conditionalPanel( condition = "output.n>=15" ,uiOutput("doseA_15")) ,
                      conditionalPanel( condition = "output.n>=16" ,uiOutput("doseA_16")) ,
                      conditionalPanel( condition = "output.n>=17" ,uiOutput("doseA_17")) ,
                      conditionalPanel( condition = "output.n>=18" ,uiOutput("doseA_18")) ,
                      conditionalPanel( condition = "output.n>=19" ,uiOutput("doseA_19")) ,
                      conditionalPanel( condition = "output.n>=20" ,uiOutput("doseA_20")) ,
                      conditionalPanel( condition = "output.n>=1", actionButton("update","Update Recommendations",icon=icon("fa fa-refresh")))
               ) , 
               column(1 , 
                      conditionalPanel( condition = "output.n>=1" ,uiOutput("doseB_1")) ,
                      conditionalPanel( condition = "output.n>=2" ,uiOutput("doseB_2")) ,
                      conditionalPanel( condition = "output.n>=3" ,uiOutput("doseB_3")) ,
                      conditionalPanel( condition = "output.n>=4" ,uiOutput("doseB_4")) ,
                      conditionalPanel( condition = "output.n>=5" ,uiOutput("doseB_5")) ,
                      conditionalPanel( condition = "output.n>=6" ,uiOutput("doseB_6")) ,
                      conditionalPanel( condition = "output.n>=7" ,uiOutput("doseB_7")) ,
                      conditionalPanel( condition = "output.n>=8" ,uiOutput("doseB_8")) ,
                      conditionalPanel( condition = "output.n>=9" ,uiOutput("doseB_9")) ,
                      conditionalPanel( condition = "output.n>=10" ,uiOutput("doseB_10")) ,
                      conditionalPanel( condition = "output.n>=11" ,uiOutput("doseB_11")) ,
                      conditionalPanel( condition = "output.n>=12" ,uiOutput("doseB_12")) ,
                      conditionalPanel( condition = "output.n>=13" ,uiOutput("doseB_13")) ,
                      conditionalPanel( condition = "output.n>=14" ,uiOutput("doseB_14")) ,
                      conditionalPanel( condition = "output.n>=15" ,uiOutput("doseB_15")) ,
                      conditionalPanel( condition = "output.n>=16" ,uiOutput("doseB_16")) ,
                      conditionalPanel( condition = "output.n>=17" ,uiOutput("doseB_17")) ,
                      conditionalPanel( condition = "output.n>=18" ,uiOutput("doseB_18")) ,
                      conditionalPanel( condition = "output.n>=19" ,uiOutput("doseB_19")) ,
                      conditionalPanel( condition = "output.n>=20" ,uiOutput("doseB_20")) 
               ) ,
               column(1 , 
                      conditionalPanel( condition = "output.n>=1" ,uiOutput("dlt1")) ,
                      conditionalPanel( condition = "output.n>=2" ,uiOutput("dlt2")) ,
                      conditionalPanel( condition = "output.n>=3" ,uiOutput("dlt3")) ,
                      conditionalPanel( condition = "output.n>=4" ,uiOutput("dlt4")) ,
                      conditionalPanel( condition = "output.n>=5" ,uiOutput("dlt5")) ,
                      conditionalPanel( condition = "output.n>=6" ,uiOutput("dlt6")) ,
                      conditionalPanel( condition = "output.n>=7" ,uiOutput("dlt7")) ,
                      conditionalPanel( condition = "output.n>=8" ,uiOutput("dlt8")) ,
                      conditionalPanel( condition = "output.n>=9" ,uiOutput("dlt9")) ,
                      conditionalPanel( condition = "output.n>=10" ,uiOutput("dlt10")) ,
                      conditionalPanel( condition = "output.n>=11" ,uiOutput("dlt11")) ,
                      conditionalPanel( condition = "output.n>=12" ,uiOutput("dlt12")) ,
                      conditionalPanel( condition = "output.n>=13" ,uiOutput("dlt13")) ,
                      conditionalPanel( condition = "output.n>=14" ,uiOutput("dlt14")) ,
                      conditionalPanel( condition = "output.n>=15" ,uiOutput("dlt15")) ,
                      conditionalPanel( condition = "output.n>=16" ,uiOutput("dlt16")) ,
                      conditionalPanel( condition = "output.n>=17" ,uiOutput("dlt17")) ,
                      conditionalPanel( condition = "output.n>=18" ,uiOutput("dlt18")) ,
                      conditionalPanel( condition = "output.n>=19" ,uiOutput("dlt19")) ,
                      conditionalPanel( condition = "output.n>=20" ,uiOutput("dlt20")) 
               )  
               
               
             ) ,
             
             br(),
             
             tableOutput("prior.ss") ,
             wellPanel(
               fluidRow(
                 column(6 , plotOutput("histplot") ) ,
                 column(6 , plotOutput("segplot") )
               )
             )
             
    ),    
    tabPanel("Simulate",
             
             
             wellPanel(fluidRow(
               column(2,offset=0, numericInput("N" , value=12 , label="Final number of patients" , min=1 , max=NA , step=1 )),
               column(2,offset=0, numericInput("S" , value=10 , label="Number of simulations" , min=1 , max=NA , step=1 )),
               column(3,offset=0, actionButton("simulate","Simulate Operating Characteristics"))
             ),
             
             fluidRow(column(4,offset=0, helpText("True p(DLT) (Drug A, Drug B)"))),
             
             fluidRow(actionButton("useprior","Use prior probabilities as truth")),
             
             fluidRow( column(1 , offset=0 ,
                              
                              conditionalPanel( condition = "input.dimA>=1 && input.dimB>=8" , numericInput("true_1_8" , value=NA , label="(1,8)" , min=0 , max=1 , step=0.05 ) ) ,
                              conditionalPanel( condition = "input.dimA>=1 && input.dimB>=7" , numericInput("true_1_7" , value=NA , label="(1,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                              conditionalPanel( condition = "input.dimA>=1 && input.dimB>=6" , numericInput("true_1_6" , value=NA , label="(1,6)" , min=0 , max=1 , step=0.05 ) ) ,
                              conditionalPanel( condition = "input.dimA>=1 && input.dimB>=5" , numericInput("true_1_5" , value=NA , label="(1,5)" , min=0 , max=1 , step=0.05 ) ) ,
                              conditionalPanel( condition = "input.dimA>=1 && input.dimB>=4" , numericInput("true_1_4" , value=NA , label="(1,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                              conditionalPanel( condition = "input.dimA>=1 && input.dimB>=3" , numericInput("true_1_3" , value=NA , label="(1,3)" , min=0 , max=1 , step=0.05 ) ) ,
                              conditionalPanel( condition = "input.dimA>=1 && input.dimB>=2" , numericInput("true_1_2" , value=NA , label="(1,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                              conditionalPanel( condition = "input.dimA>=1 && input.dimB>=1" , numericInput("true_1_1" , value=NA , label="(1,1)" , min=0 , max=1 , step=0.05 ) ) 
             ) ,
             column(1 , 
                    
                    conditionalPanel( condition = "input.dimA>=2 && input.dimB>=8" , numericInput("true_2_8" , value=NA , label="(2,8)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=2 && input.dimB>=7" , numericInput("true_2_7" , value=NA , label="(2,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=2 && input.dimB>=6" , numericInput("true_2_6" , value=NA , label="(2,6)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=2 && input.dimB>=5" , numericInput("true_2_5" , value=NA , label="(2,5)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=2 && input.dimB>=4" , numericInput("true_2_4" , value=NA , label="(2,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=2 && input.dimB>=3" , numericInput("true_2_3" , value=NA , label="(2,3)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=2 && input.dimB>=2" , numericInput("true_2_2" , value=NA , label="(2,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=2 && input.dimB>=1" , numericInput("true_2_1" , value=NA , label="(2,1)" , min=0 , max=1 , step=0.05 ) ) 
             )  ,  
             column(1 , 
                    
                    conditionalPanel( condition = "input.dimA>=3 && input.dimB>=8" , numericInput("true_3_8" , value=NA , label="(3,8)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=3 && input.dimB>=7" , numericInput("true_3_7" , value=NA , label="(3,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=3 && input.dimB>=6" , numericInput("true_3_6" , value=NA , label="(3,6)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=3 && input.dimB>=5" , numericInput("true_3_5" , value=NA , label="(3,5)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=3 && input.dimB>=4" , numericInput("true_3_4" , value=NA , label="(3,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=3 && input.dimB>=3" , numericInput("true_3_3" , value=NA , label="(3,3)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=3 && input.dimB>=2" , numericInput("true_3_2" , value=NA , label="(3,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=3 && input.dimB>=1" , numericInput("true_3_1" , value=NA , label="(3,1)" , min=0 , max=1 , step=0.05 ) ) 
             ) ,
             column(1 , 
                    
                    conditionalPanel( condition = "input.dimA>=4 && input.dimB>=8" , numericInput("true_4_8" , value=NA , label="(4,8)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=4 && input.dimB>=7" , numericInput("true_4_7" , value=NA , label="(4,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=4 && input.dimB>=6" , numericInput("true_4_6" , value=NA , label="(4,6)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=4 && input.dimB>=5" , numericInput("true_4_5" , value=NA , label="(4,5)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=4 && input.dimB>=4" , numericInput("true_4_4" , value=NA , label="(4,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=4 && input.dimB>=3" , numericInput("true_4_3" , value=NA , label="(4,3)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=4 && input.dimB>=2" , numericInput("true_4_2" , value=NA , label="(4,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=4 && input.dimB>=1" , numericInput("true_4_1" , value=NA , label="(4,1)" , min=0 , max=1 , step=0.05 ) ) 
             )  ,  
             column(1 , 
                    
                    conditionalPanel( condition = "input.dimA>=5 && input.dimB>=8" , numericInput("true_5_8" , value=NA , label="(5,8)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=5 && input.dimB>=7" , numericInput("true_5_7" , value=NA , label="(5,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=5 && input.dimB>=6" , numericInput("true_5_6" , value=NA , label="(5,6)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=5 && input.dimB>=5" , numericInput("true_5_5" , value=NA , label="(5,5)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=5 && input.dimB>=4" , numericInput("true_5_4" , value=NA , label="(5,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=5 && input.dimB>=3" , numericInput("true_5_3" , value=NA , label="(5,3)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=5 && input.dimB>=2" , numericInput("true_5_2" , value=NA , label="(5,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=5 && input.dimB>=1" , numericInput("true_5_1" , value=NA , label="(5,1)" , min=0 , max=1 , step=0.05 ) ) 
             ) ,
             column(1 , 
                    
                    conditionalPanel( condition = "input.dimA>=6 && input.dimB>=8" , numericInput("true_6_8" , value=NA , label="(6,8)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=6 && input.dimB>=7" , numericInput("true_6_7" , value=NA , label="(6,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=6 && input.dimB>=6" , numericInput("true_6_6" , value=NA , label="(6,6)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=6 && input.dimB>=5" , numericInput("true_6_5" , value=NA , label="(6,5)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=6 && input.dimB>=4" , numericInput("true_6_4" , value=NA , label="(6,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=6 && input.dimB>=3" , numericInput("true_6_3" , value=NA , label="(6,3)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=6 && input.dimB>=2" , numericInput("true_6_2" , value=NA , label="(6,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=6 && input.dimB>=1" , numericInput("true_6_1" , value=NA , label="(6,1)" , min=0 , max=1 , step=0.05 ) ) 
             )  ,  
             column(1 , 
                    conditionalPanel( condition = "input.dimA>=7 && input.dimB>=8" , numericInput("true_7_8" , value=NA , label="(7,8)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=7 && input.dimB>=7" , numericInput("true_7_7" , value=NA , label="(7,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=7 && input.dimB>=6" , numericInput("true_7_6" , value=NA , label="(7,6)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=7 && input.dimB>=5" , numericInput("true_7_5" , value=NA , label="(7,5)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=7 && input.dimB>=4" , numericInput("true_7_4" , value=NA , label="(7,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=7 && input.dimB>=3" , numericInput("true_7_3" , value=NA , label="(7,3)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=7 && input.dimB>=2" , numericInput("true_7_2" , value=NA , label="(7,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=7 && input.dimB>=1" , numericInput("true_7_1" , value=NA , label="(7,1)" , min=0 , max=1 , step=0.05 ) ) 
             ) ,
             column(1 , 
                    conditionalPanel( condition = "input.dimA>=8 && input.dimB>=8" , numericInput("true_8_8" , value=NA , label="(8,8)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=8 && input.dimB>=7" , numericInput("true_8_7" , value=NA , label="(8,7)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=8 && input.dimB>=6" , numericInput("true_8_6" , value=NA , label="(8,6)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=8 && input.dimB>=5" , numericInput("true_8_5" , value=NA , label="(8,5)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=8 && input.dimB>=4" , numericInput("true_8_4" , value=NA , label="(8,4)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=8 && input.dimB>=3" , numericInput("true_8_3" , value=NA , label="(8,3)" , min=0 , max=1 , step=0.05 ) ) ,
                    conditionalPanel( condition = "input.dimA>=8 && input.dimB>=2" , numericInput("true_8_2" , value=NA , label="(8,2)" , min=0 , max=1 , step=0.05 ) ) ,                              
                    conditionalPanel( condition = "input.dimA>=8 && input.dimB>=1" , numericInput("true_8_1" , value=NA , label="(8,1)" , min=0 , max=1 , step=0.05 ) ) 
             )    
             )) ,
             
             fluidRow( column(12,tableOutput('table'))),
             
             plotOutput("simplot")
             
             
    ) ,
    
    
    
    tabPanel("Instructions" , 
             
             helpText("Before starting, you are advised to maximise the Shiny window to the full screen width") ,   
             helpText("------------------------------------------------------------------------------") ,
             h4("Setting up the PIPE design"),
             helpText("Choose the dimensions of your dosing grid in the top-left selection boxes") ,
             helpText("Select the other design parameters") ,
             helpText("In the `Priors / Historical data' tab, choose to either enter the prior medians for the Dose Limiting Toxicities (DLTs) and the prior sample sizes or the prior numbers of DLTs / No DLTs from historical data"),
             helpText("When entering prior medians fill in the corner cells of your dosing grid first. Then either fill in the DLT prior medians for all other cells of your dosing grid OR click the `Interpolate Prior Medians' button") ,
             helpText("Note that the lowest doses are in the bottom left corner of the grid with increasing doses going up and right") ,
             fluidRow(column(12,
                             helpText("------------------------------------------------------------------------------") ,
                             h4("Conducting the trial"),
                             helpText("In the `Run Trial' tab click on the `Update Recommendations' button to get your recommended design point for the first cohort") ,
                             helpText("Complete the number of DLTs observed for this cohort and the click `Update Recommendations' again to get the recommended design points for the next cohort") ,
                             helpText("Keep on entering the number of DLTs observed at each cohort and watch your design evolve") 
             )),
             helpText("------------------------------------------------------------------------------") ,
             
             h4("Graphical output"),
             helpText("Follow the evoluation of the design in the graphs at the bottom of the page") ,
             fluidRow(column(12,helpText("The left graph shows: ")),
                      column(6, offset=2,
                             helpText("* the independent Beta posterior distributions for each combination (black lines). At the beginning of the trial these are the independent Beta priors."),
                             helpText("* the posterior distributions for each combination calculated from the monotonic contours (histograms). At the beginning of the trial these are the implied priors obtained from the monotonic contours.")
                      )),
             fluidRow(column(12,helpText("The right graph shows: ")) ,
                      
                      column(6 , offset=2 , 
                             helpText("* the recommended next design point (shaded dark grey)") ,
                             helpText("* the previous design point (shaded light grey)") ,
                             helpText("* the distribution of the maximum tolerated contour (shown by the thicknesses of the black lines)") ,                 
                             helpText("* the modal contour used by the PIPE algorithm (blue line) for the next recommended dose") ,
                             helpText("* the median contour (red line) [N.B. line is purple if modal and median contours overlap]") ,
                             helpText("* the cells that are admissible because of their proximity to the maximum tolerated contour (first box coloured green)") ,
                             helpText("* the cells that are admissible due to meeting any dose skipping constraint (second box coloured green)") ,
                             helpText("* the mean posterior probability of DLT for each combination, based on averaging over monotonic contours (p(DLT))"),
                             helpText("* the posterior probability of the DLT rate being less than the target level, based on the monotonic contours (p<theta)")
                      ))
    ) , 
    tabPanel("Background" , 
             helpText("PIPE (Product of Independent beta Probabilities dose Escalation) is a design for dual-agent phase I trials developed by Adrian Mander and Michael Sweeting at the MRC Biostatistics Unit Hub for Trials Methodology Research and the University of Cambridge in the UK."),
             helpText("Further developments of the design have been conducted in collaboration with Chris Harbron at Roche.") ,
             helpText("The Shiny implementaion of pipe.design has been jointly developed by Chris Harbron and Michael Sweeting."),
             helpText("A detailed desription of the algorithm is published in Statistics In Medicine: Volume 34, Issue 8, pages 1261-1276, 15 April 2015")
             
             
    )
  )
)
)
