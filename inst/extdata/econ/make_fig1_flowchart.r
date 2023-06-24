library(dplyr)
library(rsvg)
library(DiagrammeRsvg)

graph <- DiagrammeR::grViz("               
  # All instructions are within a large character string
digraph surveillance_diagram {    
  # 'digraph' means 'directional graph', then the graph name 
  
  # graph statement
  #################
  graph [layout = dot,
         rankdir = TB,            # layout top-to-bottom
         fontsize = 10]
  

  # nodes (circles)
  #################
  node [shape = circle,           # shape = circle
       fixedsize = true
       width = 1.3]                      
  
  Primary   [label = 'Biological\nModel'] 
  Secondary [label = 'Fishery\nModel'] 
  Tertiary  [label = 'Econometric\nModel'] 
  SC        [label = 'Estimating\nModel',
             fontcolor = darkgreen] 
  
  # edges
  #######
  Primary   -> Secondary [label = '  True\n Abundance\n (yearly)',
                          fontcolor = red,
                          color = red]
  Secondary -> Tertiary [label = '  Generated fishery\n data by year',
                          fontcolor = red,
                          color = red]
  Tertiary -> SC [label = 'Corrected\n   Abundance Index     ',
                          fontcolor = red,
                          color = red]
  Primary -> SC [label = ' Other SS3 data     ',
                          fontcolor = red,
                          color = red]                    
}
")

DPI <- 300
WidthCM <- 9

export_svg(graph) %>%
charToRaw %>%
rsvg_png(paste0(getwd(), "/inst/extdata/econ/fig1.png"),
  width = WidthCM * (DPI/2.54))
