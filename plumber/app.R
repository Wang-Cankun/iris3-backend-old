library(plumber)
#setwd('C:/Users/flyku/Documents/GitHub/iris3-backend/plumber')
setwd('/var/www/nodejs/iris3-backend/plumber')
plumber::pr("iris3.R") %>% plumber::pr_run(host="0.0.0.0",port=8000)

