# class06: R function
Siyu Xie (A16438448)

- [A first silly function](#a-first-silly-function)
- [2nd example grade() function](#2nd-example-grade-function)

All functions in R have at least 3 things: - a **name** (you get to pick
this) - input **argument** (there can be only one or loads - again your
call) - the **body** (where the work gets done, this code between the
curly brackets {})

## A first silly function

lets write a function to add some numbers. We can call it`add()`

``` r
add <- function (x) {
  y <- 10
  x+y
  
}  
```

Can I just use my new function?

``` r
add(1)
```

    [1] 11

Let’s make it a bit more flexible.

``` r
add <- function (x,y=1) {
  x+y
}  

add(x=10, y=10)
```

    [1] 20

``` r
add(10)
```

    [1] 11

## 2nd example grade() function

Write a function to grade student work.

We will start with a simple version of the problem and the following
example student vectors:

``` r
# Example input vectors to start with
student1 <- c(100, 100, 100, 100, 100, 100, 100, 90)
student2 <- c(100, NA, 90, 90, 90, 90, 97, 80)
student3 <- c(90, NA, NA, NA, NA, NA, NA, NA)
```

Start with student1

``` r
mean(student2, na.rm=TRUE)
```

    [1] 91

``` r
mean(student3, na.rm=TRUE)
```

    [1] 90

``` r
student1
```

    [1] 100 100 100 100 100 100 100  90

google told me about min() and max().

``` r
min(student1)
```

    [1] 90

``` r
which.min(student1)
```

    [1] 8

Indicate the location/index where the min occurs

``` r
student1[which.min(student1)]
```

    [1] 90

``` r
student1[-8]
```

    [1] 100 100 100 100 100 100 100

our first working snippet that drops the lowest score and calculates the
mean

``` r
mean(student1[-which.min(student1)])
```

    [1] 100

``` r
x <- student2
mean(x[-which.min(x)], na.rm=T)
```

    [1] 92.83333

``` r
x <- student3
mean(x[-which.min(x)], na.rm=T)
```

    [1] NaN

Our approach to the NA problem(missing homeworks): we can replace all NA
values with zero.

1st task is find the NA values (i.e. where are they in the vector)

``` r
x <- student2
x
```

    [1] 100  NA  90  90  90  90  97  80

``` r
is.na(x)
```

    [1] FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE

I have found the NA(TRUE) values from `is.na()` now I want to make them
equal to zero. (mask them)

``` r
x
```

    [1] 100  NA  90  90  90  90  97  80

``` r
x[is.na(x)] <- 0
```

``` r
y <- 1:5
y
```

    [1] 1 2 3 4 5

``` r
y[y>3] <-0
y
```

    [1] 1 2 3 0 0

I want to combine the `is.na(x)` with making these elements equal to
zero. And then take this “masked” (vector of student scores with NA
values as zero) and drop the lowest and get the mean.

``` r
x <- student3
x[is.na(x)] <- 0
mean(x[-which.min(x)])
```

    [1] 12.85714

Now I can turn my most awesome snippet into my first function

``` r
grade <- function(x) {
  #Make NA (missing work) equal to 0
  x[is.na(x)] <- 0
  #Drop lowest score and get the mean
  mean(x[-which.min(x)])
}
```

``` r
grade(student1)
```

    [1] 100

\#Q1. Write a function grade() to determine an overall grade from a
vector of student homework assignment scores dropping the lowest single
score. If a student misses a homework (i.e. has an NA value) this can be
used as a score to be potentially dropped. Your final function should be
adquately explained with code comments and be able to work on an example
class gradebook such as this one in CSV format:
“https://tinyurl.com/gradeinput” \[3pts\]

``` r
url <- "https://tinyurl.com/gradeinput"
  gradebook <- read.csv(url,row.names = 1)
  head(gradebook)
```

              hw1 hw2 hw3 hw4 hw5
    student-1 100  73 100  88  79
    student-2  85  64  78  89  78
    student-3  83  69  77 100  77
    student-4  88  NA  73 100  76
    student-5  88 100  75  86  79
    student-6  89  78 100  89  77

The `apply()` function in R is super useful but can be a little
confusing to begin with. Lets have a look how it works.

``` r
ans <- apply(gradebook, 1, grade)
ans
```

     student-1  student-2  student-3  student-4  student-5  student-6  student-7 
         91.75      82.50      84.25      84.25      88.25      89.00      94.00 
     student-8  student-9 student-10 student-11 student-12 student-13 student-14 
         93.75      87.75      79.00      86.00      91.75      92.25      87.75 
    student-15 student-16 student-17 student-18 student-19 student-20 
         78.75      89.50      88.00      94.50      82.75      82.75 

\#Q2. Using your grade() function and the supplied gradebook, Who is the
top scoring student overall in the gradebook?

``` r
  which.max(ans)
```

    student-18 
            18 

``` r
  max(ans)
```

    [1] 94.5

\#Q3. From your analysis of the gradebook, which homework was toughest
on students (i.e. obtained the lowest scores overall?

``` r
lowestscores <- apply(gradebook, 2, mean, na.rm=T)
lowestscores
```

         hw1      hw2      hw3      hw4      hw5 
    89.00000 80.88889 80.80000 89.63158 83.42105 

``` r
which.min(lowestscores)
```

    hw3 
      3 

``` r
min(lowestscores)
```

    [1] 80.8

\#Q4. Optional Extension: From your analysis of the gradebook, which
homework was most predictive of overall score (i.e. highest correlation
with average grade score)?

``` r
#ans
cor(gradebook$hw1,ans)
```

    [1] 0.4250204

\#mask gradebook to assign NA to 0 scores

``` r
mask <- gradebook
mask[is.na(mask)] <- 0
mask
```

               hw1 hw2 hw3 hw4 hw5
    student-1  100  73 100  88  79
    student-2   85  64  78  89  78
    student-3   83  69  77 100  77
    student-4   88   0  73 100  76
    student-5   88 100  75  86  79
    student-6   89  78 100  89  77
    student-7   89 100  74  87 100
    student-8   89 100  76  86 100
    student-9   86 100  77  88  77
    student-10  89  72  79   0  76
    student-11  82  66  78  84 100
    student-12 100  70  75  92 100
    student-13  89 100  76 100  80
    student-14  85 100  77  89  76
    student-15  85  65  76  89   0
    student-16  92 100  74  89  77
    student-17  88  63 100  86  78
    student-18  91   0 100  87 100
    student-19  91  68  75  86  79
    student-20  91  68  76  88  76

``` r
#ans
cor(mask$hw2,ans)
```

    [1] 0.176778

``` r
highestcor <- apply(mask, 2, cor, y=ans)
highestcor
```

          hw1       hw2       hw3       hw4       hw5 
    0.4250204 0.1767780 0.3042561 0.3810884 0.6325982 

``` r
which.max(highestcor)
```

    hw5 
      5 
