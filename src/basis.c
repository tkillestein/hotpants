/**
 * @file basis.c
 * @brief Kernel basis registry and dispatch mechanism.
 *
 * Manages the active kernel basis and provides initialization/cleanup.
 */

#include <stdio.h>
#include "defaults.h"
#include "globals.h"
#include "basis.h"

/* =====================================================================
   BASIS REGISTRY
   ===================================================================== */

/** Global active basis pointer */
kernel_basis_t* active_basis = NULL;

/**
 * @brief Get basis implementation for given type.
 *
 * @param[in] basis_type BASIS_TYPE_GAUSSIAN (0) or BASIS_TYPE_DELTA (1)
 * @return Pointer to basis implementation, or NULL if invalid
 */
kernel_basis_t* get_basis_for_type(int basis_type) {
  switch (basis_type) {
    case 0:  /* BASIS_TYPE_GAUSSIAN */
      return &gaussian_basis;
    case 1:  /* BASIS_TYPE_DELTA */
      return &delta_basis;
    default:
      LOG_ERROR("Unknown basis type: %d", basis_type);
      return NULL;
  }
}

/**
 * @brief Set the active basis and initialize it.
 *
 * Calls old basis cleanup() if active, then calls new basis init().
 *
 * @param[in] basis_type BASIS_TYPE_GAUSSIAN (0) or BASIS_TYPE_DELTA (1)
 * @return 0 on success, -1 on error
 */
int set_active_basis(int basis_type) {
  kernel_basis_t* new_basis;

  /* Get new basis implementation */
  new_basis = get_basis_for_type(basis_type);
  if (!new_basis) {
    LOG_ERROR("Invalid basis type: %d", basis_type);
    return -1;
  }

  /* Clean up old basis if active */
  if (active_basis) {
    if (active_basis->cleanup()) {
      LOG_ERROR("Failed to cleanup basis: %s", active_basis->name);
      return -1;
    }
  }

  /* Initialize new basis */
  if (new_basis->init() < 0) {
    LOG_ERROR("Failed to initialize basis: %s", new_basis->name);
    active_basis = NULL;
    return -1;
  }

  /* Update global nbasis for the new basis */
  new_basis->nbasis = nCompKer;
  active_basis = new_basis;

  LOG_PROGRESS("Active basis set to: %s (%d functions)", active_basis->name,
               active_basis->nbasis);

  return 0;
}

/**
 * @brief Clean up the active basis.
 *
 * @return 0 on success
 */
int cleanup_active_basis(void) {
  if (!active_basis) {
    return 0;
  }

  int ret = active_basis->cleanup();
  active_basis = NULL;
  return ret;
}
