#read xpt files, such as data from NHANES documentation

library(haven)
library(tidyverse)

setwd("S://Ye/Klotho/klotho level in serum survey/")

data <- read_xpt("SSKL_I.xpt")
write.csv(data, file = "SSKL_I.csv", row.names = FALSE)

# data$SEQN %>% length()
# data$SEQN %>% unique() %>% length()

# 2. Open a graphics device (e.g., PNG file)
pdf("HIST NHANES SSKL_I.pdf", width = 6, height = 6)

# 3. Create the histogram with custom titles
hist(data$SSKLOTH,
     main = "serum Klohto level in healthy donors (40-79 year old) \n(NHANES SSKL_I)", # Main title
     xlab = "serum Klohto (pg/mL)",             # X-axis label
     ylab = "Frequency",                 # Y-axis label
     col = "lightblue",                  # Change bar color
     border = "white",                    # Change bar border color
     breaks = 100,
)

data_median <- round(median(data$SSKLOTH, trim = 0.1),0)
data_median
data_mean <- round(mean(data$SSKLOTH, trim = 0.1),0) 
data_mean
data_n <- length(data$SSKLOTH)
data_n

# Add text annotation using the text() function
text(x = 3000, y = 200, labels = paste0("median: ", data_median, "\nmean: ", data_mean,"\nN= ", data_n), col = "black")

# Add a subtitle using the mtext() function
#mtext("This is a subtitle", side = 1, line = 4, at = 4)


# 4. Close the graphics device to save the file
dev.off()


