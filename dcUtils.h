//
//  dcUtils.h
//  cs175-final
//
//  Created by Gabe Montague on 12/16/17.
//  Copyright © 2017 cs175. All rights reserved.
//

#ifndef dcUtils_h
#define dcUtils_h

#include <stdlib.h>


namespace DC {

  double randZeroToOne() {
    return rand() / (RAND_MAX + 1.);
  }
}

#endif /* dcUtils_h */
