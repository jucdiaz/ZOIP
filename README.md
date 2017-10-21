# ZOIP
If you want to install it, you need to copy the next code into your R console

```s
if (!require("devtools")) install.packages("devtools")
if (!require("rmutil")) install.packages("rmutil")
if (!require("boot")) install.packages("boot")
if (!require("numDeriv")) install.packages("numDeriv")
if (!require("GHQp")) install.packages("GHQp")
devtools::install_github("jucdiaz/ZOIP", force = TRUE)
library(ZOIP) # Carga el paquete
```
