/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   Constants using across the program

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Constants{
    constexpr int NUM_FAM_COL = 6;
    constexpr int NUM_BIM_COL = 6;
    constexpr int NUM_MARKER_READ = 120;

    constexpr double SMALL_EPSILON = 0x1p-44; // 2^-44
}

#endif
