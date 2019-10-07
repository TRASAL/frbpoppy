c---------------------------------------------------
      real function psr_dist(x1, y1, z1, x2, y2, z2)
c---------------------------------------------------
c
c     distance between two points (Pythag C4bc)
c
      real x1,y1,x2,y2,z1,z2
      psr_dist = sqrt((((x2 - x1) ** 2) + ((y2 - y1) ** 2)) +  
     &                ((z2 - z1) **2))
      end


