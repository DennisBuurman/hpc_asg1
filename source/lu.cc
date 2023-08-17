#include <chrono>

#include <cstdio>

#include <math.h>

#include "matrix.h"


/* Global variables holding the matrix data. To complete this assignment
 * you are requested to only use arrays and access these arrays with
 * subscripts. Do not use pointers.
 */

const int max_n_elements = 131072;
const int max_n_rows = 16384;
const int test_vector_count = 5;

static double values[max_n_elements];

static int col_ind[max_n_elements];
static int row_ptr_begin[max_n_rows];
static int row_ptr_end[max_n_rows];

void print_crs(int nnz, int n_rows, int n_cols)
{
  printf("Values:\n");
  for (int i = 0; i < nnz; i++)
  {
    printf(" %.2f ", values[i]);
  }
  printf("\n");
  printf("Col_ind:\n");
  for (int i = 0; i < nnz; i++)
  {
    printf(" %d ", col_ind[i]);
  }
  printf("\n");
  printf("rb:\n");
  for (int i = 0; i < n_rows; i++)
  {
    printf(" %d ", row_ptr_begin[i]);
  }
  printf("\n");
  printf("re:\n");
  for (int i = 0; i < n_rows; i++)
  {
    printf(" %d ", row_ptr_end[i]);
  }
  printf("\n\n");
}


void print_matrix(int n_rows, int n_cols)
{
  double temp[n_cols];
  
  printf("Matrix: \n");
  
  for (int i = 0; i < n_rows; i++)
  {
    for (int i = 0; i < n_cols; i++)
    {
      temp[i] = 0.0;
    }
    for (int k = row_ptr_begin[i]; k <= row_ptr_end[i]; k++)
    {
      temp[col_ind[k]] = values[k];
    }
    for (int i = 0; i < n_cols; i++)
    {
      printf(" %.2f ", temp[i]);
    }
    printf("\n");
  }
  
  printf("\n");

}

// Refactors arrays to create new space
int collect_garbage(int &nnz, int n_rows)
{
  double new_values[max_n_elements] = {};
  
  int new_col_ind[max_n_elements] = {};
  int new_row_ptr_begin[max_n_rows] = {};
  int new_row_ptr_end[max_n_rows] = {};
  int count = 0;
  int old = nnz;
  
  // Reorder arrays in local temporary arrays
  for (int i = 0; i < n_rows; i++)
  {
    new_row_ptr_begin[i] = count;
    for (int j = row_ptr_begin[i]; j <= row_ptr_end[i]; j++)
    {
      new_values[count] = values[j];
      new_col_ind[count] = col_ind[j];
      count++;
    }
    new_row_ptr_end[i] = count - 1;
  }
  
  // Copy new arrays into global arrays
  for (int i = 0; i < max_n_elements; i++)
  {
    values[i] = new_values[i];
    col_ind[i] = new_col_ind[i];
  }
  
  for (int i = 0; i < n_rows; i++)
  {
    row_ptr_begin[i] = new_row_ptr_begin[i];
    row_ptr_end[i] = new_row_ptr_end[i];
  }
  
  nnz = count;
  
  return old-count; // Returns number of freed elements
}

// Returns the index of the pivot element in values[] and col_ind[]
int find_pivot_row(int col_index, int n_rows)
{
  int pivot_row_index = -1;
  double temp = -1.0;
  
  // Loop through all rows starting from col_index
  for (int i = col_index; i < n_rows; i++)
  {
    // Loop through row values
    for (int k = row_ptr_begin[i]; k <= row_ptr_end[i]; k++)
    {
      // Check if value in correct column
      if (col_ind[k] == col_index)
      {
        // Check if current highest found value
        if (fabs(values[k]) > temp)
        {
          temp = fabs(values[k]);
          pivot_row_index = i;
        }
      }
    }
  }
  
  return pivot_row_index;
}

// Switch row a with row b
void permutate(int row_a, int row_b)
{
  int temp_begin, temp_end;
  
  // Save positions of row a (in values[], col_ind[])
  temp_begin = row_ptr_begin[row_a];
  temp_end = row_ptr_end[row_a];
  
  // Set indexes of row b into row a
  row_ptr_begin[row_a] = row_ptr_begin[row_b];
  row_ptr_end[row_a] = row_ptr_end[row_b];
  
  // Set indexes of old row a into row b
  row_ptr_begin[row_b] = temp_begin;
  row_ptr_end[row_b] = temp_end;
}

// Eliminate values under (col_index, col_index)
void eliminate(int diagonal_index, int n_rows, int n_cols, int &nnz)
{
  double pivot = 0.0, multiplier = 0.0;
  double expanded_row[n_cols];
  double temp_row[n_cols]; // For intermediate results during dot product
  
  bool found_element;
  
  int row_index = diagonal_index;
  int temp_index; // Used in value update to index values/col_ind arrays
  int mask[n_cols];
  
  // Set the pivot
  for (int i = row_ptr_begin[diagonal_index]; i <= row_ptr_end[diagonal_index]; i++)
  {
    if (col_ind[i] == diagonal_index)
    {
      pivot = values[i];
      break;
    }
  }
  
  while(row_index < n_rows - 1)
  {
    found_element = false;
    
    // Scan rows for pivot element and set multiplier
    while (!found_element && row_index < n_rows - 1)
    {
      row_index++;
      for (int i = row_ptr_begin[row_index]; i <= row_ptr_end[row_index]; i++)
      {
        if (col_ind[i] == diagonal_index)
        {
          multiplier = values[i]/pivot;
          found_element = true;
        }
      }
    }
    
    if (found_element)
    {
      // Initialize expanded_row and temp_row
      for (int i = 0; i < n_cols; i++)
      {
        expanded_row[i] = 0.0;
        temp_row[i] = 0.0;
        mask[i] = 0.0;
      }
      
      for (int i = row_ptr_begin[row_index]; i <= row_ptr_end[row_index]; i++)
      {
        expanded_row[col_ind[i]] = values[i];
      }
      
      // Perform elimination on row
      for (int i = 0; i <= row_index; i++)
      {
        for (int j = row_ptr_begin[i]; j <= row_ptr_end[i]; j++)
        {
          // Multiply with multiplier
          if (i == diagonal_index && col_ind[j] >= i)
          {
            temp_row[col_ind[j]] -= multiplier * values[j];
          }
          else if (i == row_index)
          {
            temp_row[col_ind[j]] += values[j];
          }
        }
      }
      
      // Save multiplier in place of eliminated element
      temp_row[diagonal_index] = multiplier;
      
      int old_size = row_ptr_end[row_index] - row_ptr_begin[row_index] + 1;
      int new_size = 0;
      
      for (int i = 0; i < n_cols; i++)
      {
        if (temp_row[i] != 0.0)
        {
          new_size++;
          mask[i] = 1;
        }
      }
      
      temp_index = row_ptr_begin[row_index];
      
      if (new_size == old_size)
      {
        for (int i = 0; i < n_cols; i++)
        {
          if (mask[i] == 1)
          {
            values[temp_index] = temp_row[i];
            col_ind[temp_index] = i;
            temp_index++;
          }
        }
      }
      else if (new_size > old_size)
      {
        if (nnz + new_size < max_n_elements)
        {
          // Mark old values as garbage
          for (int i = row_ptr_begin[row_index]; i <= row_ptr_end[row_index]; i++)
          {
            values[i] = 0.0;
          }
          
          // Append adapted row at end of array
          row_ptr_begin[row_index] = nnz;
          row_ptr_end[row_index] = nnz + new_size - 1;
          temp_index = nnz;
          
          // Update values
          for (int i = 0; i < n_cols; i++)
          {
            if (mask[i] == 1)
            {
              values[temp_index] = temp_row[i];
              col_ind[temp_index] = i;
              temp_index++;
            }
          }
          
          // Increase array length
          nnz += new_size;
        }
        else
        {
          collect_garbage(nnz, n_rows);
          
          // same as if case above
          if (nnz + new_size < max_n_elements)
          {
            // Mark old values as garbage
            for (int i = row_ptr_begin[row_index]; i <= row_ptr_end[row_index]; i++)
            {
              values[i] = 0.0;
            }
            
            // Append adapted row at end of array
            row_ptr_begin[row_index] = nnz;
            row_ptr_end[row_index] = nnz + new_size - 1;
            temp_index = nnz;
            
            // Update values
            for (int i = 0; i < n_cols; i++)
            {
              if (mask[i] == 1)
              {
                values[temp_index] = temp_row[i];
                col_ind[temp_index] = i;
                temp_index++;
              }
            }
            
            // Increase array length
            nnz += new_size;
          }
        }
      }
      else // new_size < old_size
      {
        // Mark deleted values as garbage
        for (int i = row_ptr_begin[row_index] + new_size-1; i <= row_ptr_end[row_index]; i++)
        {
          values[i] = 0.0;
        }
        
        // Shrink the row_ptr range
        row_ptr_end[row_index] = row_ptr_begin[row_index] + new_size - 1;
        
        // Update values
        for (int i = 0; i < n_cols; i++)
        {
          if (mask[i] == 1)
          {
            values[temp_index] = temp_row[i];
            col_ind[temp_index] = i;
            temp_index++;
          }
        }
      }
    }
  }
}

int
main(int argc, char **argv)
{
  if (argc != 2)
    {
      fprintf(stderr, "usage: %s <filename>\n", argv[0]);
      return -1;
    }

  int nnz, n_rows, n_cols;
  bool ok(false);

  ok = load_matrix_market(argv[1], max_n_elements, max_n_rows,
                          nnz, n_rows, n_cols,
                          values, col_ind, row_ptr_begin, row_ptr_end);
  if (!ok)
    {
      fprintf(stderr, "failed to load matrix.\n");
      return -1;
    }

  /* For debugging, can be removed when implementation is finished. */
  //dump_nonzeros(n_rows, values, col_ind, row_ptr_begin, row_ptr_end);

  auto factorization_start_time = std::chrono::high_resolution_clock::now();

  /* Perform LU factorization here */
  int j = -1;
  
  for (int i = 0; i < n_rows-1; i++)
  {
    j = find_pivot_row(i, n_rows);
    if (j == -1)
    {
      continue;
    }
    
    permutate(i, j);
    eliminate(i, n_rows, n_cols, nnz);
  }

  auto factorization_end_time = std::chrono::high_resolution_clock::now();
  
  // Construct 5 b vectors
  double b_1[n_rows] = {0.0};
  double b_2[n_rows] = {0.0};
  double b_3[n_rows] = {0.0};
  double b_4[n_rows] = {0.0};
  double b_5[n_rows] = {0.0};
  
  // Construct 5 x vectors
  double x_1[n_rows];
  double x_2[n_rows];
  double x_3[n_rows];
  double x_4[n_rows];
  double x_5[n_rows];
  
  bool positive = true;
  
  for (int i = 0; i < n_rows; i++)
  {
    x_1[i] = 1.0;
    x_2[i] = 0.1;
    
    if (positive)
    {
      x_3[i] = 1.0;
      x_4[i] = 5.0;
      x_5[i] = 100.0;
    }
    else
    {
      x_3[i] = -1.0;
      x_4[i] = -5.0;
      x_5[i] = -100.0;
    }
    
    positive = !positive;
  }
  
  auto solve_start_time = std::chrono::high_resolution_clock::now();
  
  /* Compute all 5 solution vectors here */
  // Ux = c
  for (int i = 0; i < n_rows; i++)
  {
    for (int j = row_ptr_begin[i]; j <= row_ptr_end[i]; j++)
    {
      if (col_ind[j] >= i )
      {
        b_1[i] += values[j] * x_1[col_ind[j]];
        b_2[i] += values[j] * x_2[col_ind[j]];
        b_3[i] += values[j] * x_3[col_ind[j]];
        b_4[i] += values[j] * x_4[col_ind[j]];
        b_5[i] += values[j] * x_5[col_ind[j]];
      }
    }
  }
  
  // Lc = b
  for (int i = 0; i < n_rows; i++)
  {
    for (int j = row_ptr_begin[i]; j <= row_ptr_end[i]; j++)
    {
      if (col_ind[j] < i )
      {
        b_1[i] += values[j] * x_1[col_ind[j]];
        b_2[i] += values[j] * x_2[col_ind[j]];
        b_3[i] += values[j] * x_3[col_ind[j]];
        b_4[i] += values[j] * x_4[col_ind[j]];
        b_5[i] += values[j] * x_5[col_ind[j]];
      }
      else if (col_ind[j] == i)
      {
        b_1[i] += x_1[col_ind[j]];
        b_2[i] += x_2[col_ind[j]];
        b_3[i] += x_3[col_ind[j]];
        b_4[i] += x_4[col_ind[j]];
        b_5[i] += x_5[col_ind[j]];
      }
    }
  }

  auto solve_end_time = std::chrono::high_resolution_clock::now();
  
  
  double relative_errors[test_vector_count] = {0};
  
  /* Compute relative errors here */
  // TODO
  
  
  std::chrono::duration<double> factorization_elapsed_time = factorization_end_time - factorization_start_time;
  std::chrono::duration<double> solve_elapsed_time = solve_end_time - solve_start_time;
  
  
  /* Print results */
  fprintf(stdout, "%.20f\n", factorization_elapsed_time.count());
  fprintf(stdout, "%.20f\n", solve_elapsed_time.count());
  for (size_t vector_idx = 0; vector_idx < test_vector_count; ++vector_idx)
    {
      fprintf(stdout, "%.20f\n", relative_errors[vector_idx]);
    }

  return 0;
}
