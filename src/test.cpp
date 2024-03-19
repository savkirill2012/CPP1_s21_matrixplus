#include <gtest/gtest.h>

#include "S21Matrix.h"

TEST(Matrix, Constructors_1) {
  S21Matrix test;

  ASSERT_EQ(test.GetRows(), 1);
  ASSERT_EQ(test.GetCols(), 1);
  ASSERT_EQ(test(0, 0), 0);
}

TEST(Matrix, Constructors_2) {
  S21Matrix test3(5, 6);

  ASSERT_EQ(test3.GetRows(), 5);
  ASSERT_EQ(test3.GetCols(), 6);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      ASSERT_EQ(test3(i, j), 0);
    }
  }
}

TEST(Matrix, Constructors_3) {
  S21Matrix test(5, 6);
  S21Matrix test2 = test;

  ASSERT_EQ(test2.GetRows(), 5);
  ASSERT_EQ(test2.GetCols(), 6);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      ASSERT_EQ(test2(i, j), 0);
    }
  }
}

TEST(Matrix, Constructors_4) {
  S21Matrix test(5, 6);
  S21Matrix test2(std::move(test));

  ASSERT_EQ(test2.GetRows(), 5);
  ASSERT_EQ(test2.GetCols(), 6);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      ASSERT_EQ(test2(i, j), 0);
    }
  }
}

TEST(Matrix, FunctionSum) {
  S21Matrix test3(5, 6);
  S21Matrix test5(5, 6);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      test5(i, j) = j;
    }
  }

  test3.SumMatrix(test5);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      ASSERT_EQ(test3(i, j), j);
    }
  }
}

TEST(Matrix, FunctionEqual) {
  S21Matrix test3(5, 6);
  S21Matrix test5(5, 6);

  ASSERT_TRUE(test3.EqMatrix(test5));
}

TEST(Matrix, FunctionSub) {
  S21Matrix test3(5, 6);
  S21Matrix test5(5, 6);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      test5(i, j) = j;
      test3(i, j) = 1;
    }
  }

  test5.SubMatrix(test3);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      ASSERT_EQ(test5(i, j), j - 1);
    }
  }
}

TEST(Matrix, FunctionMulNum) {
  S21Matrix test3(5, 6);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      test3(i, j) = 1;
    }
  }

  test3.MulNumber(3);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      ASSERT_EQ(test3(i, j), 3);
    }
  }
}

TEST(Matrix, FunctionMulMatrix) {
  S21Matrix test3(5, 6);
  S21Matrix test5(6, 2);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      test3(i, j) = 2;
    }
  }

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 2; ++j) {
      test5(i, j) = 3;
    }
  }

  test3.MulMatrix(test5);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_EQ(test3(i, j), 36);
    }
  }
}

TEST(Matrix, FunctionTranspose) {
  S21Matrix test3(5, 6);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      test3(i, j) = 1;
    }
  }

  S21Matrix test = test3.Transpose();

  ASSERT_EQ(test.GetRows(), 6);
  ASSERT_EQ(test.GetCols(), 5);
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 5; ++j) {
      ASSERT_EQ(test(i, j), 1);
    }
  }
}

TEST(Matrix, FunctionDeterminant) {
  S21Matrix test3(3, 3);

  test3(0, 0) = 2;
  test3(0, 1) = 9;
  test3(0, 2) = 2;
  test3(1, 0) = 2;
  test3(1, 1) = 3;
  test3(1, 2) = 7;
  test3(2, 0) = 12;
  test3(2, 1) = 2;
  test3(2, 2) = 5;

  double temp = test3.Determinant();

  ASSERT_EQ((int)temp, 604);
}

TEST(Matrix, FunctionCalcComp) {
  S21Matrix test3(3, 3);

  test3(0, 0) = 2;
  test3(0, 1) = 3;
  test3(0, 2) = 8;
  test3(1, 0) = 1;
  test3(1, 1) = 9;
  test3(1, 2) = 4;
  test3(2, 0) = 13;
  test3(2, 1) = 3;
  test3(2, 2) = 5;

  S21Matrix test(test3.CalcComplements());

  ASSERT_EQ((int)test(0, 0), 33);
  ASSERT_EQ((int)test(0, 1), 47);
  ASSERT_EQ((int)test(0, 2), -114);
  ASSERT_EQ((int)test(1, 0), 9);
  ASSERT_EQ((int)test(1, 1), -94);
  ASSERT_EQ((int)test(1, 2), 33);
  ASSERT_EQ((int)test(2, 0), -60);
  ASSERT_EQ((int)test(2, 1), 0);
  ASSERT_EQ((int)test(2, 2), 15);
}

TEST(Matrix, FunctionInversivMatrix) {
  S21Matrix test3(3, 3);

  test3(0, 0) = 2;
  test3(0, 1) = 3;
  test3(0, 2) = 1;
  test3(1, 0) = 0;
  test3(1, 1) = 5;
  test3(1, 2) = 2;
  test3(2, 0) = 3;
  test3(2, 1) = 1;
  test3(2, 2) = 0;

  S21Matrix test = test3.InverseMatrix();

  ASSERT_EQ((int)test(0, 0), 2);
  ASSERT_EQ((int)test(0, 1), -1);
  ASSERT_EQ((int)test(0, 2), -1);
  ASSERT_EQ((int)test(1, 0), -6);
  ASSERT_EQ((int)test(1, 1), 3);
  ASSERT_EQ((int)test(1, 2), 4);
  ASSERT_EQ((int)test(2, 0), 15);
  ASSERT_EQ((int)test(2, 1), -7);
  ASSERT_EQ((int)test(2, 2), -10);
}

TEST(Matrix, Operations) {
  S21Matrix test3(5, 6);
  S21Matrix test5(5, 6);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      test3(i, j) = 2;
      test5(i, j) = 3;
    }
  }

  test3 += test5;
  test3 = test3 + test3;

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_EQ(test3(i, j), 10);
    }
  }

  test3 *= 2;

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_EQ(test3(i, j), 20);
    }
  }

  test3 = test3 - test5;
  if (test3 == test3) {
    test3 -= test5;
  }

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_EQ(test3(i, j), 14);
    }
  }

  test3 = test5;

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 2; ++j) {
      ASSERT_EQ(test3(i, j), 3);
    }
  }
}

TEST(Matrix, FunctionSetters) {
  S21Matrix test3(5, 6);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 6; ++j) {
      test3(i, j) = 1;
    }
  }

  test3.SetRows(3);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 5; ++j) {
      ASSERT_EQ(test3(i, j), 1);
    }
  }

  test3.SetCols(3);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      ASSERT_EQ(test3(i, j), 1);
    }
  }

  test3.SetRows(5);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (i > 2) {
        ASSERT_EQ(test3(i, j), 0);
      } else {
        ASSERT_EQ(test3(i, j), 1);
      }
    }
  }

  test3.SetCols(5);

  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      if (i > 2 || j > 2) {
        ASSERT_EQ(test3(i, j), 0);
      } else {
        ASSERT_EQ(test3(i, j), 1);
      }
    }
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}