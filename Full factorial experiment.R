Sys.setlocale("LC_ALL","Russian_Russia")

get.averages <- function(minValues, maxValues)
{
	avg <- (minValues + maxValues) / 2
	
	cat("Нулевые значения факторов:\n")
	for (i in 1:length(avg))
	{
		cat("x", i, " = ", avg[i], "  ", sep = "")
	}
	cat("\n\n")
	
	return(avg)
}

get.deviation <- function(minValues, maxValues)
{
	dev <- (maxValues - minValues) / 2
	
	cat("Отклонения факторов:\n")
	for (i in 1:length(dev))
	{
		cat("Δx", i, " = ", dev[i], "  ", sep = "")
	}
	cat("\n\n")

	return(dev)
}

new.expand.grid <- function(vec, nrep)
{
	do.call(expand.grid, rep(list(vec), nrep))
}

get.planning.matrix <- function(count)
{
	x <- c("-1", "1")
	matrix <- new.expand.grid(x, count)
	for (i in 1:count)
	{
		name <- paste("x", as.character(i), sep = "")
		names(matrix)[i]<- name
	}
	
	x0 <- rep(1, nrow(matrix))
	data <- data.frame(x0, matrix)
	
	for (i in seq(from = 2, to = count))
	{
		bases <- seq(1:count)
		while (TRUE)
		{
			name <- ""
			newCol <- rep(1, nrow(matrix))
      		for (j in seq(from = 1, to = i))
			{
				col <- as.numeric(as.character(unlist(as.vector(data[bases[j]+1]), use.names=FALSE)))
				newCol <- newCol * col
				name <- paste(name, paste("x", as.character(bases[j]), sep = ""), sep = "")
			}

			data <- data.frame(data, newCol)
			names(data)[length(data)]<- name
		
			if (bases[1] + (i-1) == count)
			{
				break
			}

			for (j in seq(from = i, to = 1, by = -1))
			{
				if (bases[j] - (j - 1) < count - (i - 1))
				{
					bases[j] <- bases[j] + 1
					if (bases[j] - j < count - i)
					{
						for(k in seq(from = j, to = i))
						{
							bases[k+1] = bases[k] + 1
						}
					}
				
					break
				}
			}
		}
	}
	
	cat("Матрица планирования:\n")
	print(data)
	cat("\n")
	
	return(data)
}

get.average.results <- function(y)
{
	avg <- 0
	for(i in 1:length(y))
	{
		avg <- avg + y[i]
	}
	avg <- unlist(avg / length(y), recursive = TRUE, use.names = FALSE)
	
	return(avg)
}

get.dispersion <- function(y, avg_y)
{
	dispersion <- rep(0, length(avg_y))
	for (i in 1:length(avg_y))
	{
		for (j in 1:length(y))
		{
			dispersion[i] <- dispersion[i] + ((y[i,j] - avg_y[i]) ** 2)
		}
		dispersion[i] <- dispersion[i] / (length(y) - 1)
	}
	
	return(dispersion)
}

show.dispersions.matrix <- function(y, avg_y, dispersion)
{
	matrix <- data.frame(y, avg_y)
	names(matrix)[length(matrix)]<- "y"
	matrix <- data.frame(matrix, dispersion)
	names(matrix)[length(matrix)]<- "S"
	
	cat("Средние значения и выборочные дисперсии результатов опытов:\n")
	print(matrix)
	cat("\n")
}

get.Cochrans.criterion <- function(y, avg_y, dispersion)
{
	maxDispersion <- max(dispersion)
	sumDispersion <- sum(dispersion)
	criterion <- maxDispersion / sumDispersion
	
	cat("Критерий Кохрена: ")
	cat(criterion, "\n")
	
	return(criterion)
}

check.Cochrans.criterion <- function(criterion, y)
{
	f <- length(y) - 1
	N <- as.character(nrow(y))
	
	table = read.table("Cochrans criterion.txt", header=TRUE)
	tableValue <- table[N,f]
	
	cat("Табличное значение критерия Кохрена: ")
	cat(tableValue, "\n")
	
	if (criterion < tableValue)
	{
		cat("Эксперимент воспроизводим\n\n")
		return(TRUE)
	}
	else
	{
		cat("Эксперимент невоспроизводим\n\n")
		return(FALSE)
	}
}

show.source.regression.equation <- function(x, count)
{
	cat("Исходное уравнение регрессии:\n")
	cat("y = b0")
	coefs_names <- c("b0")

	for (i in 1:count)
	{
		cat(" + b", i, "*", "x", i, sep = "")
		coefs_names <- append(coefs_names, paste("b", as.character(i), sep = ""))
	}

	
	for (i in seq(from = 2, to = count))
	{
		bases <- seq(1:count)
		while (TRUE)
		{
			cat(" + ")
			
			coef <- "b"
			name <- ""
      		for (j in seq(from = 1, to = i))
			{
				coef <- paste(coef, as.character(bases[j]), sep = "")
				name <- paste(name, "*", paste("x", as.character(bases[j]), sep = ""), sep = "")
			}
			
			cat(paste(coef, name, sep = ""))
			coefs_names <- append(coefs_names, coef)
		
			if (bases[1] + (i-1) == count)
			{
				break
			}

			for (j in seq(from = i, to = 1, by = -1))
			{
				if (bases[j] - (j - 1) < count - (i - 1))
				{
					bases[j] <- bases[j] + 1
					if (bases[j] - j < count - i)
					{
						for(k in seq(from = j, to = i))
						{
							bases[k+1] = bases[k] + 1
						}
					}
				
					break
				}
			}
		}
	}
	
	cat("\n\n")

	return(coefs_names)
}

get.coefficients  <- function(x, avg_y, coefs_names)
{
	coefs <- solve(matrix(as.numeric(as.matrix(x)), ncol = ncol(x)), avg_y)
	
	cat("Коэффициенты уравнения регрессии:\n")
	for (i in 1:length(coefs))
	{
		cat(coefs_names[i], " = ", coefs[i], "\n", sep = "")
	}
	cat("\n")

	return(coefs)
}

get.reproducibility.dispersion <- function(dispersion)
{
	avg_dispersion <- mean(dispersion)
	
	cat("Дисперсия воспроизводимости: ")
	cat(avg_dispersion, "\n\n")
	
	return(avg_dispersion)
}

get.Students.criterion <- function(coefs, avg_dispersion, y, coefs_names)
{
	criteria <- rep(1, length(coefs))
	numerator <- nrow(y) * length(y)
	for (i in 1:length(coefs))
	{
		criteria[i] <- criteria[i] * abs(coefs[i]) * sqrt(numerator / avg_dispersion)
	}
	
	cat("Критерий Стьюдента:\n")
	for (i in 1:length(criteria))
	{
		cat("t(", coefs_names[i], ") = ", criteria[i], "\n", sep = "")
	}
	cat("\n")
	
	return(criteria)
}

check.Students.criterion <- function(criteria, y, coefs, coefs_names)
{
	f <- nrow(y) * (length(y) - 1)
	
	table = read.table("Students criterion.txt", header=FALSE)
	tableValue <- table[f, 1]
	
	cat("Табличное значение критерия Стьюдента: ")
	cat(tableValue, "\n")
	
	cat("Значимые коэффициенты:\n")
	
	significance <- rep(0, length(criteria))
	for (i in 1:length(criteria))
	{
		if (criteria[i] > tableValue)
		{
			cat(coefs_names[i], " = ", coefs[i], "\n", sep = "")
			
			significance[i] <- 1
		}
	}
	cat("\n")
	
	return(significance)
}

get.significance.coefficients <- function(coefs, significance)
{
	significance.coefs <- c()
	for (i in 1:length(coefs))
	{
		if (significance[i] == 1)
		{
			significance.coefs <- append(significance.coefs, coefs[i])
		}
	}
	
	return(significance.coefs)
}

show.regression.equation <- function(coefs, significance, x)
{
	equation <- "y = "
	significance_coefs <- c()
	signs <- c()
	for (i in 1:length(coefs))
	{
		if (significance[i] == 1)
		{
			if (i == 1)
			{
				name <- as.character(abs(coefs[i]))
			}
			else
			{
				name <- as.character(abs(coefs[i]))
				name <- paste(name, colnames(x)[i], sep = "")
			}
			significance_coefs <- append(significance_coefs, name)
			if (coefs[i] > 0)
			{
				signs <- append(signs, 1)
			}
			else
			{
				signs <- append(signs, 0)
			}
		}
	}
	
	if (length(significance_coefs) > 0)
	{
		if (signs[1] > 0)
		{
			equation <- paste(equation, significance_coefs[1], sep = "")
		}
		else
		{
			equation <- paste(equation, significance_coefs[1], sep = "- ")
		}
		
		if (length(significance_coefs) > 1)
		{
			for (i in 2:length(significance_coefs))
			{
				if (signs[i] > 0)
				{
					equation <- paste(equation, significance_coefs[i], sep = " + ")
				}
				else
				{
					equation <- paste(equation, significance_coefs[i], sep = " - ")
				}
			}
		}
	}
	else
	{
		equation <- "отсутствует"
	}
	
	cat("Полученное уравнение регрессии:\n")
	cat(equation, "\n\n")
}

get.result.by.equation <- function(significance, x, significance.coefs)
{
	significance.matrix <- matrix(as.numeric(as.matrix(x)), ncol = ncol(x))
	for (i in seq(from = ncol(x), to = 1, by = -1))
	{
		if (significance[i] == 0)
		{	
			significance.matrix <- matrix(significance.matrix[,-i], ncol = ncol(significance.matrix) - 1)
		}
	}
	
	results <- rep(0, nrow(x))
	for (i in 1:nrow(x))
	{
		for (j in 1:length(significance.coefs))
		{
			results[i] <- results[i] + significance.coefs[j] * significance.matrix[i,j]
		}
	}
	
	cat("Результаты опытов для полученного уравнения:\n")
	cat(results, "\n\n")
	
	return(results)
}

get.residual.dispersion <- function(significance.coefs, results, y, avg_y)
{
	multiplier <- length(y) / (nrow(y) - length(significance.coefs))
	
	residual.dispersion <- 0
	for (i in 1:nrow(y))
	{
		residual.dispersion <- residual.dispersion + ((results[i] - avg_y[i]) ** 2)
	}
	residual.dispersion <- residual.dispersion * multiplier
	
	cat("Остаточная дисперсия (дисперсия адекватности): ")
	cat(residual.dispersion, "\n\n")
	
	return(residual.dispersion)
}

get.Fishers.criterion <- function(avg_dispersion, residual.dispersion)
{
	criterion <- residual.dispersion / avg_dispersion
	
	cat("Критерий Фишера: ")
	cat(criterion, "\n")

	return(criterion)
}

check.Fishers.criterion <- function(criterion, significance.coefs, y)
{
	k1 <- paste("X", as.character(nrow(y) - length(significance.coefs)), sep = "")
	k2 <- as.character(nrow(y) * (length(y) - 1))
	
	table = read.table("Fishers criterion.txt", header=TRUE)
	tableValue <- table[k2, k1]
	
	cat("Табличное значение критерия Фишера: ")
	if (length(tableValue) > 0)
	{
		cat(tableValue, "\n")
	
		if (criterion < tableValue)
		{
			cat("Уравнение регрессии адекватно\n\n")
			return(TRUE)
		}
		else
		{
			cat("Уравнение регрессии неадекватно\n\n")
			return(FALSE)
		}
	}
	else
	{
		cat("отсутствует в таблице\n\n")
		return(FALSE)
	}
}

full.factorial.experiment <- function(minZ, maxZ, y)
{
	isCorrect <- FALSE
	if (is.vector(minZ) & is.vector(maxZ))
	{
		if (length(minZ) == length(maxZ))
		{
			z <- get.averages(minZ, maxZ)
			dev_z <- get.deviation(minZ, maxZ)
			
			x <- get.planning.matrix(length(z))

			if (nrow(x) == nrow(y))
			{
				avg_y <- get.average.results(y)
				dispersion <- get.dispersion(y, avg_y)
				show.dispersions.matrix(y, avg_y, dispersion)

				cochrans_criterion <- get.Cochrans.criterion(y, avg_y, dispersion)
				reproduced <- check.Cochrans.criterion(cochrans_criterion, y)
				
				if (reproduced)
				{
					coefs_names <- show.source.regression.equation(x, length(z))
					coefs <- get.coefficients(x, avg_y, coefs_names)
					
					avg_dispersion <- get.reproducibility.dispersion(dispersion)
					students_criterion <- get.Students.criterion(coefs, avg_dispersion, y, coefs_names)
					significance <- check.Students.criterion(students_criterion, y, coefs, coefs_names)
					
					significance.coefs <- get.significance.coefficients(coefs, significance)
					show.regression.equation(coefs, significance, x)
					
					if (length(significance.coefs) > 0)
					{
						results <- get.result.by.equation(significance, x, significance.coefs)
						residual.dispersion <- get.residual.dispersion(significance.coefs, results, y, avg_y)
						fishers.criterion <- get.Fishers.criterion(avg_dispersion, residual.dispersion)
						adequately <- check.Fishers.criterion(fishers.criterion, significance.coefs, y)
					}
				}
				
				
				isCorrect <- TRUE
			}
		}
	}

	if (!isCorrect)
		cat("Некорректные параметры.\n")
}


#Source data
min <- c(220, 0.8, 350)
max <- c(300, 1.6, 700)

y1 <- c(60, 250, 125, 220, 225, 235, 125, 225)
y2 <- c(120, 175, 160, 225, 195, 155, 115, 170)
y3 <- c(100, 220, 95, 205, 210, 275, 180, 250)
y4 <- c(80, 230, 80, 175, 190, 245, 145, 260)
y <- data.frame(y1, y2, y3, y4)


#Run experiment
full.factorial.experiment(min, max, y)