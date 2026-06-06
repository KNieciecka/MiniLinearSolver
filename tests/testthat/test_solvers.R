# testy jednostkowe

# rozklad lu
test_that("Poprawność rozkładu LU: PA = LU", {
  
  A <- matrix(c(2, 1, 1, 4, -6, 0, -2, 7, 2), nrow = 3)
  
  wynik <- lu_decompose(A)
  
  PA <- wynik$P %*% A
  LU <- wynik$L %*% wynik$U
  
  expect_equal(as.numeric(PA), as.numeric(LU))
})

# pivot
test_that("Poprawność zamiany wierszy w rozkładzie LU", {
  A_pivot <- matrix(c(0,1,-1,2,-2,1,1,3,2), nrow=3)
  
  wynik <- lu_decompose(A_pivot)
  
  expect_equal(as.numeric(wynik$P %*% A_pivot), as.numeric(wynik$L %*% wynik$U))
  
  expect_false(isTRUE(all.equal(wynik$P, diag(3))))
})

# macierz jednostkowa
test_that("Poprawność rozkładu Lu dla macierzy jednostkowej", {
  I <- diag(4)
  
  wynik <- lu_decompose(I)
  
  expect_equal(as.numeric(wynik$L), as.numeric(I))
  expect_equal(as.numeric(wynik$U), as.numeric(I))
  expect_equal(lu_det(I), 1)
})

# zerowa kolumna
test_that("Poprawność rozkładu LU dla macierzy z zerową kolumną", {
  A_zer <- matrix(c(0, 0, 0, 2, 4, 1, 5, 3, 2), nrow = 3, ncol = 3)
  b_test <- c(1, 2, 3)
  
  expect_error(suppressWarnings(lu_solve(A_zerowa, b_test)))
})

# wyznacznik
test_that("Poprawność obliczania wyznacznika", {
  A <- matrix(c(2, 4,-2,1,-6,7,1,0,2), nrow=3)
  
  expect_equal(lu_det(A), det(A))
  
  A_osob <- matrix(c(1, 2, 2, 4), nrow=2)
  
  expect_equal(suppressWarnings(lu_det(A_osob)), 0)
})

# rozwiązywanie LU
test_that("Poprawność rozwiązywania układu metodą LU", {
  A <- matrix(c(2, 4,-2,1,-6,7,1,0,2), nrow=3)
  b <- c(4, -2, 7)
  
  expect_equal(lu_solve(A, b), solve(A, b))
})

# układ sprzeczny
test_that("Poprawność metody Lu dla układu sprzecznego",{
  A_sprzecz <- matrix(c(1, 1, 1, 1), nrow = 2 )
  b_sprzecz <- c(2, 3)

  expect_error(suppressWarnings(lu_solve(A_sprzecz, b_sprzecz)))
})

# ukłąd nieoznaczony
test_that("Poprawność metody LU dla układu nieoznaczonego",{
  A_nieoz <- matrix(c(1, 2, 1, 2), nrow = 2)
  b_nieoz <- c(2, 4)

  expect_error(suppressWarnings(lu_solve(A_nieoz, b_nieoz)))
})

#male wartosci
test_that("Poprawność rozwiązywania układu LU dla małych wartości", {
  A_m <- matrix(c(1e-8, 3e-8,2e-8, 4e-8), nrow=2)
  b_m <- c(1e-8, 2e-8)
  
  expect_equal(lu_solve(A_m, b_m), solve(A_m, b_m))
})

# rozwiązywanie SOR
test_that("Poprawność rozwiązywania układu metodą SOR", {
  A <- matrix(c(10,-1,2,0,-1,11,-1,3,2,-1,10,-1,0,3,-1,8), nrow=4)
  b <- c(6, 25, -11, 15)
  
  wynik_sor <- sor_solve(A, b, omega = 1.1)
  wynik_r <- solve(A, b)
  
  expect_equal(as.numeric(wynik_sor$x), as.numeric(wynik_r), tolerance = 1e-6)
})

# nie zbieżne SOR
test_that("Poprawność metody SOR dla macierzy rozbieżnej", {
  A_rozb <- matrix(c(1, 10,10, 1), nrow=2)
  b <- c(1, 1)
  limit <- 100
  
  wynik_sor <- suppressWarnings(sor_solve(A_rozb, b, max_iter = limit))
  
  expect_equal(wynik_sor$iterations, limit)
})

# układ sprzeczny
test_that("Poprawność metody SOR dla układu sprzecznego", {
  A_sprzecz <- matrix(c(1, 1, 1, 1), nrow = 2)
  b_sprzecz <- c(2, 3)
  
  limit <- 100
  
  wynik_sor <- suppressWarnings(sor_solve(A_sprzecz, b_sprzecz, max_iter = limit))
  
  expect_equal(wynik_sor$iterations, limit)
})

# układ nieoznaczony
test_that("Poprawność metody SOR dla układu nieoznaczonego", {
  A_nieoz <- matrix(c(1, 2, 1, 2), nrow = 2)
  b_nieoz <- c(2, 4)
  limit <- 100
  
  wynik_sor <- suppressWarnings(sor_solve(A_nieoz, b_nieoz, max_iter = limit))
  
  lewa_str <- as.numeric(A_nieoz %*% wynik_sor$x)
  
  expect_equal(lewa_str, b_nieoz, tolerance = 1e-6)
  
  expect_true(wynik_sor$iterations < limit)
})

# errory
test_that("Poprawność wyrzucania błędów", {
  A_niekwad <- matrix(1:6, nrow=2)
  b <- c(1, 2)
  
  expect_error(lu_decompose(A_niekwad))
  expect_error(lu_solve(A_niekwad, b))
  expect_error(sor_solve(A_niekwad,b))

  expect_error(lu_det("to jest jakis tekst"))
})