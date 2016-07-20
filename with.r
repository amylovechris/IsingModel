with(mtcars, {
	nokeepstats <- summary(mpg)
	keepstats <<- summary(mpg)
})
nokeepstats
keepstats

