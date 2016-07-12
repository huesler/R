
# http://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf


################################################################################################
################################################################################################
################################################################################################

#Syntax - Helpful conventions for wrangling ----------------------------------------

dplyr::tbl_df(iris)
dplyr::glimpse(iris)
dplyr::glimpse(iris)
dplyr::%>%

x %>% f(y)  is the same as  f(x, y)
y %>% f(x, ., z)  is the same as  f(x, y, z )

iris %>%
group_by(Species) %>%
summarise(avg = mean(Sepal.Width)) %>%
arrange(avg)


#Subset Observations (Rows) --------------------------------------------------------

dplyr::filter(iris, Sepal.Length > 7)

dplyr::distinct(iris)

dplyr::sample_frac(iris, 0.5, replace = TRUE)

dplyr::sample_n(iris, 10, replace = TRUE)

dplyr::slice(iris, 10:15)

dplyr::top_n(storms, 2, date)


#Subset Variables (Columns) --------------------------------------------------------

dplyr::select(iris, Sepal.Width, Petal.Length, Species)

#Helper functions for select - ?select

select(iris, contains("."))

select(iris, ends_with("Length"))

select(iris, everything())

select(iris, matches(".t."))

select(iris, num_range("x", 1:5))

select(iris, one_of(c("Species", "Genus")))

select(iris, starts_with("Sepal"))

select(iris, Sepal.Length:Petal.Width)

select(iris, -Species)


#Reshaping Data - Change the layout of a data set ----------------------------------

tidyr::gather(cases, "year", "n", 2:4)

tidyr::separate(storms, date, c("y", "m", "d"))

tidyr::spread(pollution, size, amount)

tidyr::unite(data, col, ..., sep)



dplyr::data_frame(a = 1:3, b = 4:6)

dplyr::arrange(mtcars, mpg)

dplyr::arrange(mtcars, desc(mpg))

dplyr::rename(tb, y = year)


################################################################################################
################################################################################################
################################################################################################

#Summarise Data --------------------------------------------------------------------
dplyr::summarise(iris, avg = mean(Sepal.Length))

dplyr::summarise_each(iris, funs(mean))

dplyr::count(iris, Species, wt = Sepal.Length)


dplyr::first

dplyr::last

dplyr::nth

dplyr::n

dplyr::n_distinct


IQR

min

max

mean

median

var

sd


#Make New Variables ----------------------------------------------------------------

dplyr::mutate(iris, sepal = Sepal.Length + Sepal. Width)

dplyr::mutate_each(iris, funs(min_rank))

dplyr::transmute(iris, sepal = Sepal.Length + Sepal. Width)


dplyr::lead

dplyr::lag

dplyr::dense_rank

dplyr::min_rank

dplyr::percent_rank

dplyr::row_number

dplyr::ntile

dplyr::between

dplyr::cume_dist

dplyr::cumall

dplyr::cumany

dplyr::cummean


cumsum

cummax

cummin

cumprod

pmax

pmin


#Group Data ------------------------------------------------------------------------

dplyr::group_by(iris, Species)

dplyr::ungroup(iris)

iris %>% group_by(Species) %>% summarise(…)

iris %>% group_by(Species) %>% mutate(…)

#Combine Data Sets -----------------------------------------------------------------

#Mutating Joins
dplyr::left_join(a, b, by = "x1")

dplyr::right_join(a, b, by = "x1")

dplyr::inner_join(a, b, by = "x1")

dplyr::full_join(a, b, by = "x1")


#Filtering Joins
dplyr::semi_join(a, b, by = "x1")

dplyr::anti_join(a, b, by = "x1")


#Set Operations
dplyr::intersect(y, z)

dplyr::union(y, z)

dplyr::setdiff(y, z)


#Binding
dplyr::bind_rows(y, z)

dplyr::bind_cols(y, z)


