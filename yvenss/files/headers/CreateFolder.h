/*
 * CreateFolder.h
 *
 *  Created on: Mar 11, 2017
 *      Author: ary
 */

#ifndef HEADERS_CREATEFOLDER_H_
#define HEADERS_CREATEFOLDER_H_

#include <sys/stat.h>
#include <sys/types.h>

#include "PathConstants.h"

    class CreateFolder {
    public:
        CreateFolder(string path);
        virtual ~CreateFolder();
    };

#endif /* HEADERS_CREATEFOLDER_H_ */
