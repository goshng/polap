
################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with 
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################


g2 = g %>% flip(27,31) %>% rotate(29) +
  geom_cladelab(node=36, size=1, label="Oleaceae", align=TRUE, offset = 0.15, offset.text=0.005) +
  geom_cladelab(node=14, label="Gesneriaceae", align=TRUE, offset = 0.15, offset.text=0.005) +
  geom_cladelab(node=25, label="Plantaginaceae", align=TRUE, offset = 0.15, offset.text=0.005) +
  geom_cladelab(node=13, label="Phrymaceae", align=TRUE, offset = 0.15, offset.text=0.005) +
  geom_cladelab(node=32, label="Orobanchaceae", align=TRUE, offset = 0.15, offset.text=0.005) +
  geom_cladelab(node=27, label="Lamiaceae", align=TRUE, offset = 0.15, offset.text=0.005)

p2 %>% insert_left(g2, width=1.8) 