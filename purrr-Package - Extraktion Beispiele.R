
#Reference Manual - purrr-Package - Extraktion Beispiele:
#https://cran.r-project.org/web/packages/purrr/purrr.pdf (copy-paste)
######################################################################
######################################################################


# 1 accumulate . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 3

1:3 %>% accumulate(‘+‘)
1:10 %>% accumulate_right(‘*‘)
# From Haskell's scanl documentation
1:10 %>% accumulate(max, .init = 5)
# Simulating stochastic processes with drift
## Not run:
library(dplyr)
library(ggplot2)
rerun(5, rnorm(100)) %>%
set_names(paste0("sim", 1:5)) %>%
map(~ accumulate(., ~ .05 + .x + .y)) %>%
map_df(~ data_frame(value = .x, step = 1:100), .id = "simulation") %>%
ggplot(aes(x = step, y = value)) +
geom_line(aes(color = simulation)) +
ggtitle("Simulations of a random walk with drift")
## End(Not run)


# 2 along . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4

x <- 1:5
rep_along(x, 1:2)
rep_along(x, 1)
list_along(x)


# 3 array-coercion . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 4

# We create an array with 3 dimensions
x <- array(1:12, c(2, 2, 3))
# A full margin for such an array would be the vector 1:3. This is
# the default if you don't specify a margin
# Creating a branch along the full margin is equivalent to
# as.list(array) and produces a list of size length(x):
array_branch(x) %>% str()
# A branch along the first dimension yields a list of length 2
# with each element containing a 2x3 array:
array_branch(x, 1) %>% str()
# A branch along the first and third dimensions yields a list of
# length 2x3 whose elements contain a vector of length 2:
array_branch(x, c(1, 3)) %>% str()
# Creating a tree from the full margin creates a list of lists of
# lists:
array_tree(x) %>% str()
# The ordering and the depth of the tree are controlled by the
# margin argument:
array_tree(x, c(3, 1)) %>% str()


# 4 as_function . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 5

as_function(~ . + 1)
as_function(1)
as_function(c("a", "b", "c"))
as_function(c("a", "b", "c"), .null = NA)


# 5 as_vector . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 6

# Supply the type either with a string:
as.list(letters) %>% as_vector("character")
# Or with a vector mold:
as.list(letters) %>% as_vector(character(1))
# Vector molds are more flexible because they also specify the
# length of the concatenated vectors:
list(1:2, 3:4, 5:6) %>% as_vector(integer(2))
# Note that unlike vapply(), as_vector() never adds dimension
# attributes. So when you specify a vector mold of size > 1, you
# always get a vector and not a matrix


# 6 at_depth . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 7

l1 <- list(
obj1 = list(
prop1 = list(param1 = 1:2, param2 = 3:4),
prop2 = list(param1 = 5:6, param2 = 7:8)
),
obj2 = list(
prop1 = list(param1 = 9:10, param2 = 11:12),
prop2 = list(param1 = 13:14, param2 = 15:16)
)
)
# In the above list, "obj" is level 1, "prop" is level 2 and "param"
# is level 3. To apply sum() on all params, we map it at depth 3:
l1 %>% at_depth(3, sum)
# map() lets us pluck the elements prop1/param2 in obj1 and obj2:
l1 %>% map(c("prop1", "param2")) %>% str()
# But what if we want to pluck all param2 elements? Then we need to
# act at a lower level:
l1 %>% at_depth(2, "param2") %>% str()
# at_depth can be used in a complementary way with other purrr
# functions to make them operate at a lower level
l2 <- list(
obj1 = list(
prop1 = list(c(1, 2), c(3, 4), c(5, 6)),
prop2 = list(c("a", "b"), c("c", "d"), c("e", "f"))
),
obj2 = list(
	prop1 = list(c(10, 20), c(30, 40), c(50, 60)),
prop2 = list(c("A", "B"), c("C", "D"), c("E", "F"))
)
)
# Here we ask pmap() to map paste() simultaneously over all
# elements of the objects at the second level. paste() is thus
# effectively mapped at level 3.
l2 %>% at_depth(2, pmap, paste, sep = " / ")	


# 7 bare-type-predicates . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 9


# 8 by_row . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 10

# ..f should be able to work with a list or a data frame. As it
# happens, sum() handles data frame so the following works:
mtcars %>% by_row(sum)
# Other functions such as mean() may need to be adjusted with one
# of the lift_xy() helpers:
mtcars %>% by_row(lift_vl(mean))
# To run a function with invoke_rows(), make sure it is variadic (that
# it accepts dots) or that .f's signature is compatible with the
# column names
mtcars %>% invoke_rows(.f = sum)
mtcars %>% invoke_rows(.f = lift_vd(mean))
# invoke_rows() with cols collation is equivalent to plyr::mdply()
p <- expand.grid(mean = 1:5, sd = seq(0, 1, length = 10))
p %>% invoke_rows(.f = rnorm, n = 5, .collate = "cols")
## Not run:
p %>% plyr::mdply(rnorm, n = 5) %>% dplyr::tbl_df()
## End(Not run)
# To integrate the result as part of the data frame, use rows or
# cols collation:
mtcars[1:2] %>% by_row(function(x) 1:5)
mtcars[1:2] %>% by_row(function(x) 1:5, .collate = "rows")
mtcars[1:2] %>% by_row(function(x) 1:5, .collate = "cols")



# 9 by_slice . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 11

# Here we fit a regression model inside each slice defined by the
# unique values of the column "cyl". The fitted models are returned
# in a list-column.
mtcars %>%
slice_rows("cyl") %>%
by_slice(partial(lm, mpg ~ disp))
# by_slice() is especially useful in combination with map().
# To modify the contents of a data frame, use rows collation. Note
# that unlike dplyr, Mutating and summarising operations can be
# used indistinctly.
# Mutating operation:
df <- mtcars %>% slice_rows(c("cyl", "am"))
df %>% by_slice(dmap, ~ .x / sum(.x), .collate = "rows")
# Summarising operation:
df %>% by_slice(dmap, mean, .collate = "rows")
# Note that mapping columns within slices is best handled by dmap():
df %>% dmap(~ .x / sum(.x))
df %>% dmap(mean)
# If you don't need the slicing variables as identifiers, switch
# .labels to FALSE:
mtcars %>%
slice_rows("cyl") %>%
by_slice(partial(lm, mpg ~ disp), .labels = FALSE) %>%
flatten() %>%
map(coef)


# 10 compose . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 13

not_null <- compose(‘!‘, is.null)
not_null(4)
not_null(NULL)
add1 <- function(x) x + 1
compose(add1, add1)(8)


# 11 conditional-map . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 14

# Convert factors to characters
iris %>%
map_if(is.factor, as.character) %>%
str()
# Specify which columns to map with a numeric vector of positions:
mtcars %>% map_at(c(1, 4, 5), as.character) %>% str()
# Or with a vector of names:
mtcars %>% map_at(c("cyl", "am"), as.character) %>% str()
list(x = rbernoulli(100), y = 1:100) %>%
transpose() %>%
map_if("x", ~ update_list(., y = ~ y * 100)) %>%
transpose() %>%
simplify_all()

# 12 contains . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 15

x <- list(1:10, 5, 9.9)
x %>% contains(1:10)
x %>% contains(3)


# 13 cross_n . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 15

# We build all combinations of names, greetings and separators from our
# list of data and pass each one to paste()

data <- list(
id = c("John", "Jane"),
greeting = c("Hello.", "Bonjour."),
sep = c("! ", "... ")
)
data %>%
cross_n() %>%
map(lift(paste))

# cross_n() returns the combinations in long format: many elements,
# each representing one combination. With cross_d() we'll get a
# data frame in long format: crossing three objects produces a data
# frame of three columns with each row being a particular
# combination. This is the same format that expand.grid() returns.
args <- data %>% cross_d()
# In case you need a list in long format (and not a data frame)
# just run as.list() after cross_d()
args %>% as.list()
# This format is often less pratical for functional programming
# because applying a function to the combinations requires a loop
out <- vector("list", length = nrow(args))
for (i in seq_along(out))
out[[i]] <- map(args, i) %>% invoke(paste, .)
out
# It's easier to transpose and then use invoke_map()
args %>% transpose() %>% map_chr(~ invoke(paste, .))
# Unwanted combinations can be filtered out with a predicate function
filter <- function(x, y) x >= y
cross2(1:5, 1:5, .filter = filter) %>% str()
# To give names to the components of the combinations, we map
# setNames() on the product:
seq_len(3) %>%
cross2(., ., .filter = ‘==‘) %>%
map(setNames, c("x", "y"))
# Alternatively we can encapsulate the arguments in a named list
# before crossing to get named components:
seq_len(3) %>%
list(x = ., y = .) %>%
cross_n(.filter = ‘==‘)

# 14 detect . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 17

is_even <- function(x) x %% 2 == 0
3:10 %>% detect(is_even)
3:10 %>% detect_index(is_even)
3:10 %>% detect(is_even, .right = TRUE)
3:10 %>% detect_index(is_even, .right = TRUE)


# 15 dmap . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 18

# dmap() always returns a data frame:
dmap(mtcars, summary)
# dmap() also supports sliced data frames:
sliced_df <- mtcars[1:5] %>% slice_rows("cyl")
sliced_df %>% dmap(mean)
sliced_df %>% dmap(~ .x / max(.x))
# This is equivalent to the combination of by_slice() and dmap()
# with 'rows' collation of results:
sliced_df %>% by_slice(dmap, mean, .collate = "rows")


# 16 every . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 19

x <- list(0, 1, TRUE)
x %>% every(identity)
x %>% some(identity)
y <- list(0:10, 5.5)
y %>% every(is.numeric)
y %>% every(is.integer)


# 17 flatten . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 20

#These functions remove a level hierarchy from a list. They are similar to unlist, only ever remove
# a single layer of hierarchy, and are type-stable so you always know what the type of the output is.
x <- rerun(2, sample(4))
x
x %>% flatten()
x %>% flatten_int()
# You can use flatten in conjunction with map
x %>% map(1L) %>% flatten_int()
# But it's more efficient to use the typed map instead.
x %>% map_int(1L)





# 18 get-attr . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 21

factor(1:3) %@% "levels"
mtcars %@% "class"


# 19 head_while . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 22
pos <- function(x) x >= 0
head_while(5:-5, pos)
tail_while(5:-5, negate(pos))
big <- function(x) x > 100
head_while(0:10, big)
tail_while(0:10, big)



# 20 invoke . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 22

# Invoke a function with a list of arguments
invoke(runif, list(n = 10))
# Invoke a function with named arguments
invoke(runif, n = 10)
# Combine the two:
invoke(paste, list("01a", "01b"), sep = "-")
# That's more natural as part of a pipeline:
list("01a", "01b") %>%
invoke(paste, ., sep = ".")
# Invoke a list of functions, each with different arguments
invoke_map(list(runif, rnorm), list(list(n = 10), list(n = 5)))
# Or with the same inputs:
invoke_map(list(runif, rnorm), list(list(n = 5)))
invoke_map(list(runif, rnorm), n = 5)
# Or the same function with different inputs:
invoke_map("runif", list(list(n = 5), list(n = 10)))
# Or as a pipeline
list(m1 = mean, m2 = median) %>% invoke_map(x = rcauchy(100))
list(m1 = mean, m2 = median) %>% invoke_map_dbl(x = rcauchy(100))
# Note that you can also match by position by explicitly omitting ‘.x‘.
# This can be useful when the argument names of the functions are not
# identical
list(m1 = mean, m2 = median) %>%
invoke_map(, rcauchy(100))
# If you have pairs of function name and arguments, it's natural
# to store them in a data frame:
if (requireNamespace("dplyr", quietly = TRUE)) {
df <- dplyr::data_frame(
f = c("runif", "rpois", "rnorm"),
params = list(
list(n = 10),
list(n = 5, lambda = 10),
list(n = 10, mean = -3, sd = 10)
)
)
df
invoke_map(df$f, df$params)
}


# 21 is_empty . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 24

is_empty(NULL)
is_empty(list())
is_empty(list(NULL))


# 22 is_formula . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 25

x <- disp ~ am
is_formula(x)


# 23 keep . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 25

rep(10, 10) %>%
map(sample, 5) %>%
keep(function(x) mean(x) > 6)
# Or use a formula
rep(10, 10) %>%
map(sample, 5) %>%
keep(~ mean(.x) > 6)
# Using a string instead of a function will select all list elements
# where that subelement is TRUE
x <- rerun(5, a = rbernoulli(1), b = sample(10))
x
x %>% keep("a")
x %>% discard("a")


# 24 lift . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 26

### Lifting from ... to list(...) or c(...)
x <- list(x = c(1:100, NA, 1000), na.rm = TRUE, trim = 0.9)
lift_dl(mean)(x)
# Or in a pipe:
mean %>% lift_dl() %>% invoke(x)
# You can also use the lift() alias for this common operation:
lift(mean)(x)
# Default arguments can also be specified directly in lift_dl()
list(c(1:100, NA, 1000)) %>% lift_dl(mean, na.rm = TRUE)()
# lift_dl() and lift_ld() are inverse of each other.
# Here we transform sum() so that it takes a list
fun <- sum %>% lift_dl()
fun(list(3, NA, 4, na.rm = TRUE))
# Now we transform it back to a variadic function
fun2 <- fun %>% lift_ld()
fun2(3, NA, 4, na.rm = TRUE)
# It can sometimes be useful to make sure the lifted function's
# signature has no named parameters, as would be the case for a
# function taking only dots. The lifted function will take a list
# or vector but will not match its arguments to the names of the
# input. For instance, if you give a data frame as input to your
# lifted function, the names of the columns are probably not
# related to the function signature and should be discarded.
lifted_identical <- lift_dl(identical, .unnamed = TRUE)
mtcars[c(1, 1)] %>% lifted_identical()
mtcars[c(1, 2)] %>% lifted_identical()
#
### Lifting from c(...) to list(...) or ...
# Some functions such as mean() take an atomic vector. It is often
# useful to transform them to functions taking a list. In the
# following example, we lift mean() to apply it to each row of a
# data frame. This works because a row is essentially a list of
# length-1 vectors:
mtcars %>% by_row(lift_vl(mean))
# In other situations we need the vector-valued function to take a
# variable number of arguments as with pmap(). This is a job for
# lift_vd():
pmap(mtcars, lift_vd(mean))
# lift_vd() will collect the arguments and concatenate them to a
# vector before passing them to ..f. You can add a check to assert
# the type of vector you expect:
lift_vd(tolower, .type = character(1))("this", "is", "ok")
#
### Lifting from list(...) to c(...) or ...
# cross_n() normally takes a list of elements and returns their
# cartesian product. By lifting it you can supply the arguments as
# if it was a function taking dots:
cross <- lift_ld(cross_n)
out1 <- cross_n(list(a = 1:2, b = c("a", "b", "c")))
out2 <- cross(a = 1:2, b = c("a", "b", "c"))
identical(out1, out2)
# This kind of lifting is sometimes needed for function
# composition. An example would be to use pmap() with a function
# that takes a list. In the following, we use some() on each row of
# a data frame to check they each contain at least one element
# satisfying a condition:
mtcars %>% pmap(lift_ld(some, partial(‘<‘, 200)))
# Default arguments for ..f can be specified in the call to
# lift_ld()
lift_ld(cross_n, .filter = ‘==‘)(1:3, 1:3) %>% str()
# Here is another function taking a list and that we can update to
# take a vector:
glue <- function(l) {
if (!is.list(l)) stop("not a list")
l %>% invoke(paste, .)
}
## Not run:
letters %>% glue() # fails because glue() expects a list
## End(Not run)
letters %>% lift_lv(glue)() # succeeds

# 25 lmap . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 29

# Let's write a function that returns a larger list or an empty list
# depending on some condition. This function also uses the names
# metadata available in the attributes of the list-element
maybe_rep <- function(x) {
n <- rpois(1, 2)
out <- rep_len(x, n)
if (length(out) > 0) {
names(out) <- paste0(names(x), seq_len(n))
}
out
}
# The output size varies each time we map f()
x <- list(a = 1:4, b = letters[5:7], c = 8:9, d = letters[10])
x %>% lmap(maybe_rep)
# We can apply f() on a selected subset of x
x %>% lmap_at(c("a", "d"), maybe_rep)
# Or only where a condition is satisfied
x %>% lmap_if(is.character, maybe_rep)

# A more realistic example would be a function that takes discrete
# variables in a dataset and turns them into disjunctive tables, a
# form that is amenable to fitting some types of models.
# A disjunctive table contains only 0 and 1 but has as many columns
# as unique values in the original variable. Ideally, we want to
# combine the names of each level with the name of the discrete
# variable in order to identify them. Given these requirements, it
# makes sense to have a function that takes a data frame of size 1
# and returns a data frame of variable size.
disjoin <- function(x, sep = "_") {
name <- names(x)
x <- as.factor(x[[1]])
out <- lapply(levels(x), function(level) {
as.numeric(x == level)
})
names(out) <- paste(name, levels(x), sep = sep)
dplyr::as_data_frame(out)
}
# Now, we are ready to map disjoin() on each categorical variable of a
# data frame:
iris %>% lmap_if(is.factor, disjoin)
mtcars %>% lmap_at(c("cyl", "vs", "am"), disjoin)

# 26 map . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 31

1:10 %>%
map(rnorm, n = 10) %>%
map_dbl(mean)
# Or use an anonymous function
1:10 %>%
map(function(x) rnorm(10, x))
# Or a formula
1:10 %>%
map(~ rnorm(10, .x))
# A more realistic example: split a data frame into pieces, fit a
# model to each piece, summarise and extract R^2
mtcars %>%
split(.$cyl) %>%
map(~ lm(mpg ~ wt, data = .x)) %>%
map(summary) %>%
map_dbl("r.squared")
# Use map_lgl(), map_dbl(), etc to reduce to a vector.
# * list
mtcars %>% map(sum)
# * vector
mtcars %>% map_dbl(sum)
# If each element of the output is a data frame, use
# map_df to row-bind them together:
mtcars %>%
split(.$cyl) %>%
map(~ lm(mpg ~ wt, data = .x)) %>%
map_df(~ as.data.frame(t(as.matrix(coef(.)))))
# (if you also want to preserve the variable names see
# the broom package)


# 27 map2 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 33
################################################# Auch unter map2 zu finden!!!
#map2(.x, .y, .f, ...)
#map2_lgl(.x, .y, .f, ...)
#map2_int(.x, .y, .f, ...)
#map2_dbl(.x, .y, .f, ...)
#map2_chr(.x, .y, .f, ...)
#map2_df(.x, .y, .f, ..., .id = NULL)
#pmap(.l, .f, ...)
#pmap_lgl(.l, .f, ...)
#pmap_int(.l, .f, ...)
#pmap_dbl(.l, .f, ...)
#pmap_chr(.l, .f, ...)
#pmap_df(.l, .f, ..., .id = NULL)
#walk2(.x, .y, .f, ...)
#pwalk(.l, .f, ...)

x <- list(1, 10, 100)
y <- list(1, 2, 3)
map2(x, y, ~ .x + .y)
# Or just
map2(x, y, ‘+‘)
# Split into pieces, fit model to each piece, then predict
by_cyl <- mtcars %>% split(.$cyl)
mods <- by_cyl %>% map(~ lm(mpg ~ wt, data = .))
map2(mods, by_cyl, predict)


# 28 negate . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 35

x <- transpose(list(x = 1:10, y = rbernoulli(10)))
x %>% keep("y") %>% length()
x %>% keep(negate("y")) %>% length()
# Same as
x %>% discard("y") %>% length()


# 29 null-default . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 36

1 %||% 2
NULL %||% 2


# 30 partial . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 36

# Partial is designed to replace the use of anonymous functions for
# filling in function arguments. Instead of:
compact1 <- function(x) discard(x, is.null)
# we can write:
compact2 <- partial(discard, .p = is.null)
# and the generated source code is very similar to what we made by hand
compact1
compact2
# Note that the evaluation occurs "lazily" so that arguments will be
# repeatedly evaluated
f <- partial(runif, n = rpois(1, 5))
f
f()
f()
# You can override this by saying .lazy = FALSE
f <- partial(runif, n = rpois(1, 5), .lazy = FALSE)
f
f()
f()
# This also means that partial works fine with functions that do
# non-standard evaluation
my_long_variable <- 1:10
plot2 <- partial(plot, my_long_variable)
plot2()
plot2(runif(10), type = "l")


# 31 prepend . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 37

x <- as.list(1:3)
x %>% append("a")
x %>% prepend("a")
x %>% prepend(list("a", "b"), before = 3)


# 32 rbernoulli . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 38

rbernoulli(10)
rbernoulli(100, 0.1)


# 33 rdunif . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 39

table(rdunif(1e3, 10))


# 34 reduce . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 39

1:3 %>% reduce(‘+‘)
1:10 %>% reduce(‘*‘)
5 %>%
replicate(sample(10, 5), simplify = FALSE) %>%
reduce(intersect)
x <- list(c(0, 1), c(2, 3), c(4, 5))
x %>% reduce(c)
x %>% reduce_right(c)
# Equivalent to:
x %>% rev() %>% reduce(c)
# Use init when you want reduce to return a consistent type when
# given an empty lists
list() %>% reduce(‘+‘)
list() %>% reduce(‘+‘, .init = 0)


# 35 rerun . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 40

10 %>% rerun(rnorm(5))
10 %>%
rerun(x = rnorm(5), y = rnorm(5)) %>%
map_dbl(~ cor(.x$x, .x$y))


# 36 safely . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 41

safe_log <- safely(log)
safe_log(10)
safe_log("a")
list("a", 10, 100) %>%
map(safe_log) %>%
transpose()
# This is a bit easier to work with if you supply a default value
# of the same type and use the simplify argument to transpose():
safe_log <- safely(log, otherwise = NA_real_)
list("a", 10, 100) %>%
map(safe_log) %>%
transpose() %>%
simplify_all()
# To replace errors with a default value, use possibly().
list("a", 10, 100) %>%
map_dbl(possibly(log, NA_real_))


# 37 scalar-type-predicates . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 42

#These predicates check for a given type and whether the vector is "scalar", that is, of length 1.


# 38 set_names . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 43

set_names(1:4, c("a", "b", "c", "d"))
# If the second argument is ommitted a vector is named with itself
set_names(letters[1:5])


# 39 slice_rows . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 43

#slice_rows() is equivalent to dplyr’s group_by() command but it takes a vector of column names
#or positions instead of capturing column names with special evaluation. unslice() removes the
#slicing attributes.



# 40 splice . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 44

inputs <- list(arg1 = "a", arg2 = "b")
# splice() concatenates the elements of inputs with arg3
splice(inputs, arg3 = c("c1", "c2")) %>% str()
list(inputs, arg3 = c("c1", "c2")) %>% str()
c(inputs, arg3 = c("c1", "c2")) %>% str()


# 41 split_by . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 44

#Auch hier zu finden
#split_by(.x, .f, ...)
#order_by(.x, .f, ...)
#sort_by(.x, .f, ...)

l1 <- transpose(list(x = sample(10), y = 1:10))
l1
l1 %>% order_by("x")
l1 %>% sort_by("x")
l2 <- rerun(5, g = sample(2, 1), y = rdunif(5, 10))
l2 %>% split_by("g") %>% str()
l2 %>% split_by("g") %>% map(. %>% map("y"))

# 42 transpose . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 45

x <- rerun(5, x = runif(1), y = runif(5))
x %>% str()
x %>% transpose() %>% str()
# Back to where we started
x %>% transpose() %>% transpose() %>% str()
# transpose() is useful in conjunction with safely() & quietly()
x <- list("a", 1, 2)
y <- x %>% map(safely(log))
y %>% str()
y %>% transpose() %>% str()
# Use simplify_all() to reduce to atomic vectors where possible
x <- list(list(a = 1, b = 2), list(a = 3, b = 4), list(a = 5, b = 6))
x %>% transpose()
x %>% transpose() %>% simplify_all()


# 43 type-predicates . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 46

is_list(x)
is_atomic(x)
is_vector(x)
is_numeric(x)
is_integer(x)
is_double(x)
is_character(x)
is_logical(x)
is_null(x)
is_function(x)


# 44 update_list . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 47

x <- list(x = 1:10, y = 4)
update_list(x, z = 10)
update_list(x, z = ~ x + y)


# 45 when . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 48

1:10 %>%
when(
sum(.) <= 50 ~ sum(.),
sum(.) <= 100 ~ sum(.)/2,
~ 0
)
1:10 %>%
when(
sum(.) <= x ~ sum(.),
sum(.) <= 2*x ~ sum(.)/2,
~ 0,
x = 60
)
iris %>%
subset(Sepal.Length > 10) %>%
when(
nrow(.) > 0 ~ .,
~ iris %>% head(10)
)
iris %>%
head %>%
when(nrow(.) < 10 ~ .,
~ stop("Expected fewer than 10 rows."))


