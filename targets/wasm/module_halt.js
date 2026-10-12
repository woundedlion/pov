/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/* global Module */

(function () {
  const onAbort = Module['onAbort'];
  Module['onAbort'] = function (reason) {
    Module['HS_MODULE_DEAD'] = true;
    onAbort?.call(Module, reason);
  };
})();
