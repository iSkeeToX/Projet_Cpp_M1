#include <iostream>
#include <vector>

void dfs(int i, int j, std::vector<std::vector<int>>& grid, std::vector<std::vector<bool>>& visited, int currentLabel) {
    if (i < 0 || i >= grid.size() || j < 0 || j >= grid[0].size() || visited[i][j] || grid[i][j] == 0) {
        return;
    }

    visited[i][j] = true;
    grid[i][j] = currentLabel;

    // Explore neighbors
    dfs(i + 1, j, grid, visited, currentLabel);
    dfs(i - 1, j, grid, visited, currentLabel);
    dfs(i, j + 1, grid, visited, currentLabel);
    dfs(i, j - 1, grid, visited, currentLabel);
}

void findConnectedComponents(std::vector<std::vector<int>>& grid) {
    int rows = grid.size();
    int cols = grid[0].size();
    std::vector<std::vector<bool>> visited(rows, std::vector<bool>(cols, false));
    int currentLabel = 1;  // Start labeling from 1

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (!visited[i][j] && grid[i][j] == 1) {
                dfs(i, j, grid, visited, currentLabel);
                currentLabel++;
            }
        }
    }
}

int main() {
    // Example usage
    std::vector<std::vector<int>> grid = {
        {1, 0, 1, 0},
        {1, 1, 0, 1},
        {0, 1, 1, 0},
        {0, 0, 0, 1}
    };

    findConnectedComponents(grid);

    // Display the labeled grid
    for (const auto& row : grid) {
        for (int cell : row) {
            std::cout << cell << " ";
        }
        std::cout << "\n";
    }

    return 0;
}

