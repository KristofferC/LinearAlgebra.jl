#!/usr/bin/env bash

while read -r OS ARCH ROOTFS_IMAGE_NAME ROOTFS_TAG ROOTFS_ARCH ROOTFS_HASH; do
    # Skip empty lines and comments
    [[ -z "${OS}" || "${OS}" == \#* ]] && continue
      export OS ARCH ROOTFS_IMAGE_NAME ROOTFS_TAG ROOTFS_ARCH ROOTFS_HASH
      buildkite-agent pipeline upload ./.buildkite/runtests.yml
done < ".buildkite/runtests-matrix.txt"
