#include "s21_matrix.h"

// calloc

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int error = OK;
  if (rows > 0 && columns > 0) {
    if (result) {
      result->rows = rows;
      result->columns = columns;
      result->matrix = calloc(rows, sizeof(double *));
      for (int i = 0; i < rows; i++) {
        result->matrix[i] = (double *)calloc(columns, sizeof(double));
      }
      if (!result->matrix) {
        free(result->matrix);
        error = ERROR;
      }
    } else {
      error = ERROR;
    }
  } else {
    error = ERROR;
  }
  return error;
}

void s21_remove_matrix(matrix_t *A) {
  if (s21_check_matrix(A) == 0) {
    for (int i = 0; i < A->rows; i++) {
      free(A->matrix[i]);
    }
    free(A->matrix);
    A->matrix = NULL, A->rows = 0, A->columns = 0;
  }
}

// enum function_result_eq {SUCCESS = 1, FAILURE = 0}

// FIXED

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int error = SUCCESS;
  if (s21_check_matrix(A) == ERROR || s21_check_matrix(B) == ERROR) {
    error = FAILURE;
  } else if (A->rows == B->rows && A->columns == B->columns) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) {
          error = FAILURE;
          break;
        }
      }
      if (error == FAILURE) break;
    }
  } else {
    error = FAILURE;
  }
  return error;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) == OK && s21_check_matrix(B) == OK &&
      s21_check_size(A, B) == OK) {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else {
      error = ERROR;
    }
  } else if (s21_check_matrix(A) == ERROR && s21_check_matrix(B) == ERROR) {
    error = ERROR;
  } else {
    error = ERROR_CALC;
  }
  return error;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) == OK && s21_check_matrix(B) == OK &&
      s21_check_size(A, B) == OK) {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    } else {
      error = ERROR;
    }
  } else if (s21_check_matrix(A) == ERROR && s21_check_matrix(B) == ERROR) {
    error = ERROR;
  } else {
    error = ERROR_CALC;
  }
  return error;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) == OK) {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] * number;
        }
      }
    } else {
      error = ERROR;
    }
  } else if (s21_check_matrix(A) == ERROR) {
    error = ERROR;
  } else {
    error = ERROR_CALC;
  }
  return error;
}

// MEMORY FAIL RESULT FAIL

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) == OK && s21_check_matrix(B) == OK &&
      A->columns == B->rows) {
    if (s21_create_matrix(A->rows, B->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          for (int k = 0; k < A->columns; k++) {
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
          }
        }
      }
    } else {
      error = ERROR;
    }
  } else if (s21_check_matrix(A) != OK || s21_check_matrix(B) != OK) {
    error = ERROR;
  } else {
    error = ERROR_CALC;
  }
  return error;
}

// FIXED

int s21_transpose(matrix_t *A, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) == OK) {
    s21_create_matrix(A->columns, A->rows, result);
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  } else {
    error = ERROR;
  }
  return error;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int error = OK;
  if ((s21_check_matrix(A) != OK)) {
    error = ERROR;
  } else if (A->rows != A->columns) {
    error = ERROR_CALC;
  } else {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      if (A->rows == 1) {
        result->matrix[0][0] = A->matrix[0][0];
      } else {
        double det_tmp = 0;
        for (int i = 0; i < A->rows; i++) {
          for (int j = 0; j < A->columns; j++) {
            matrix_t tmp;
            s21_create_matrix(A->rows - 1, A->columns - 1, &tmp);
            s21_tmp_matrix(A, &tmp, i, j);
            s21_determinant(&tmp, &det_tmp);
            result->matrix[i][j] = det_tmp * pow(-1, i + j);
            s21_remove_matrix(&tmp);
            det_tmp = 0;
          }
        }
      }
    }
  }
  return error;
}

int s21_determinant(matrix_t *A, double *result) {
  int error = OK;
  double det_tmp = 0;
  *result = 0;
  if (s21_check_matrix(A) != OK) {
    error = ERROR;
  } else if (A->rows != A->columns) {
    error = ERROR_CALC;
  } else {
    if (A->rows == 1) {
      *result = A->matrix[0][0];
    } else if (A->rows == 2) {
      *result =
          A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
    } else {
      for (int i = 0; i < A->columns; i++) {
        matrix_t det;
        s21_create_matrix(A->rows - 1, A->columns - 1, &det);
        s21_tmp_matrix(A, &det, 0, i);
        s21_determinant(&det, &det_tmp);
        *result += A->matrix[0][i] * det_tmp * pow(-1, i);
        s21_remove_matrix(&det);
      }
    }
  }
  return error;
}

// MEMORY FAIL RESULT FAIL

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int error = OK;
  if (s21_check_matrix(A) != OK) {
    error = ERROR;
  } else if (A->rows != A->columns) {
    error = ERROR_CALC;
  } else {
    double res;
    if (s21_determinant(A, &res) == OK) {
      if (fabs(res) == 0) {
        error = ERROR_CALC;
      } else {
        matrix_t det, inverse;
        s21_create_matrix(A->rows, A->columns, result);
        s21_calc_complements(A, &det);
        s21_transpose(&det, &inverse);
        s21_mult_number(&inverse, 1 / res, result);
        s21_remove_matrix(&det);
        s21_remove_matrix(&inverse);
      }
    } else {
      error = ERROR_CALC;
    }
  }
  return error;
}

// extra function

void s21_tmp_matrix(matrix_t *A, matrix_t *tmp, int r, int c) {
  for (int i = 0, I_tmp = 0; i < A->rows; i++) {
    if (i != r) {
      for (int j = 0, J_tmp = 0; j < A->columns; j++) {
        if (j != c) {
          tmp->matrix[I_tmp][J_tmp] = A->matrix[i][j];
          J_tmp++;
        }
      }
      I_tmp++;
    }
  }
}

int s21_check_matrix(matrix_t *A) {
  int error = OK;
  if (A == NULL || A->matrix == NULL || A->rows < 1 || A->columns < 1 ||
      !A->matrix) {
    error = ERROR;
  }
  return error;
}

int s21_check_size(matrix_t *A, matrix_t *B) {
  int error = ERROR;
  if (A->rows == B->rows && A->columns == B->columns) {
    error = OK;
  }
  return error;
}
