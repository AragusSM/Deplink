library(plotly)

#Generate a plotky with massage
DepLink_Plotly_Message <- function(selected_type = c("drug", "gene"), position = c("below", "above")){
  
msg <- paste0("Please select one ", selected_type," from the table ",position)
height <- 1
width <- 1
data <- data.frame(msg, height,width)

fig <- plot_ly(data, x = ~width, y = ~height, type = 'scatter',
               mode = 'text', text = ~msg, textposition = 'middle',
               textfont = list(color = '#000000', size = 16))
fig <- fig %>% layout(title = '',
                      xaxis = list(showticklabels = FALSE,title = "",zeroline = F,showgrid = F),
                      yaxis = list(showticklabels = FALSE,title = "", zeroline = F,showgrid = F))

return(fig)
}


#DepLink_Plotly_Message(selected_type = "gene",position = "above")
