#include <iostream>
#include <iomanip>
#include <vector>
#include <optional>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <utility>

using Matrix = std::vector<std::vector<int64_t>>;
using OptionalMatrix = std::optional<Matrix>;
using Point = std::pair<int64_t, int64_t>;
using Cover = std::tuple<Point, Point, Point>;
using CoverVector = std::vector<Cover>;
using OptionalCoverVector = std::optional<CoverVector>;

std::ostream& operator << (std::ostream& out, const Point& point)
{
  out << "(" << point.first << ", " << point.second << ")";
  return out;
}

std::ostream& operator << (std::ostream& out, const Cover& cover)
{
  out << "<";
  out << std::get<0>(cover) << ", ";
  out << std::get<1>(cover) << ", ";
  out << std::get<2>(cover);
  out << ">";

  return out;
}

void printUsage()
{
  std::cout << "Usage: ./mn_table <row-count> <coulmn-cout>\n";
}
int64_t digitCount(int64_t number)
{
  return std::log2(number) + 1;
}

void printMatrix(const Matrix& matrix)
{
  const int64_t width = matrix.size();
  const int64_t height = matrix[0].size();
  const int64_t cellLen = digitCount(width * height);

  for(int64_t i = 0; i < width; i++) {
    for(int64_t j = 0; j < height; j++) {
      std::cout << std::setw(cellLen) << matrix[i][j];
    }
    std::cout << "\n";
  }
}

int64_t edgeCount(int64_t n, int64_t m) 
{
  if(n == 0 || m == 0) {
    return 0;
  }

  return (n - 1) * m + n * (m - 1);
}

OptionalMatrix constructTable(int64_t n, int64_t m)
{
  if(n == 0 || m == 0) {
    return std::nullopt;
  }

  Matrix table;
  table.resize(n);
  for(int64_t i = 0; i < n; i++) {
    table[i].resize(m);
  }

  for(int64_t i = 0; i < n; i++) {
    for(int64_t j = i; j < m; j++) {
      table[i][j] = edgeCount(i, j);
    }
  }

  return table;
}

Cover returnCover(Point point)
{
  Point p1 {point.first + 1, point.second};
  Point p2 {point.first, point.second};
  Point p3 {point.first, point.second + 1};

  return std::make_tuple(std::move(p1), std::move(p2), std::move(p3));
}

void coverRightmostColumnFromOutside(CoverVector& resultCover, int64_t nCopy, int64_t mCopy)
{
  for(int64_t i = mCopy; i > 1; i -= 2) {
    resultCover.emplace_back(Cover{{nCopy, i}, {nCopy, i - 1}, {nCopy, i - 2}});
  }
}

void coverLowerRowFromOutside(CoverVector& resultCover, int64_t, int64_t mCopy)
{  
  for(int64_t i = 0; i < mCopy - 1; i += 2) {
    resultCover.emplace_back(Cover{{i, 0}, {i + 1, 0}, {i + 2, 0}});
  }
}

void coverRightmostBottomCell(CoverVector& resultCover, int64_t, int64_t mCopy)
{
  resultCover.emplace_back(Cover{{mCopy - 1, 0}, {mCopy, 0}, {mCopy, + 1}});
}

OptionalCoverVector doTwoEdgeCover(int64_t n, int64_t m)
{
  // -----------------------------------
  // To do edge cover, the number of edges should be divisible by 2
  const auto edgeCountLocal = edgeCount(n, m);
  std::cerr << "\nedgeCount(" << n << ", " << m << ") = " << edgeCountLocal << "\n";
  // -----------------------------------

  // -----------------------------------
  bool isGraphAllowed = !(edgeCountLocal % 2);
  if(!isGraphAllowed) {
    std::cout << "Number of edges isn't divisible by 2\nCan't perform two-edge cover.\n";
    return std::nullopt;
  }
  // -----------------------------------

  // -----------------------------------
  // Determine, if there's exists redundunt column for outer cover
  bool isRedundantColumn = !(n % 2) && !(m % 2);
  std::cout << std::boolalpha << "Is redundunt column exists: " << isRedundantColumn << std::endl;
  // -----------------------------------

  // -----------------------------------
  // Make internal cover
  // For each point return possible cover from it
  auto resultCover = CoverVector{};
  for(int64_t i = 0; i <= n - 1; i++) {
    for(int64_t j = 0; j <= m - 1; j++) {
      resultCover.emplace_back(returnCover({i, j}));
    }
  }
  // -----------------------------------

  // -----------------------------------
  int64_t nCopy = n;
  int64_t mCopy = m;
  if(isRedundantColumn) {
    nCopy -= 1;
    mCopy -= 1;
  }
  // -----------------------------------

  // -----------------------------------
  // Cover rightmost column from outside
  coverRightmostColumnFromOutside(resultCover, nCopy, mCopy);
  // -----------------------------------

  // -----------------------------------
  // Cover lower row from outside
  coverLowerRowFromOutside(resultCover, nCopy,  mCopy);
  // -----------------------------------

  // -----------------------------------
  // Cover rightmost bottom cell
  coverRightmostBottomCell(resultCover, nCopy, mCopy);
  // -----------------------------------
  
  return resultCover;
}

int main(int argc, char *argv[])
{
  if(argc != 3) {
    printUsage();
    return 1;
  }
  
  int64_t n = std::stoul(argv[1]);
  int64_t m = std::stoul(argv[2]);;

  auto optCoverVector = doTwoEdgeCover(n, m);
  auto coverVector = CoverVector{};
  if(optCoverVector.has_value()) {
    coverVector = optCoverVector.value();
  } else {
    std::cout << "Cover for " << m << "x" << n << " table doesn't exists\n";
  }

  for(const auto& cover : coverVector) {
    std::cout << cover << std::endl;
  }

  return 0;
}
