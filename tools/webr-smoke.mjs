import { createServer } from "node:http";
import { readFile } from "node:fs/promises";
import { extname, resolve, sep } from "node:path";
import { WebR } from "webr";

const repositoryRoot = resolve("_site");
const mimeTypes = new Map([
  [".gz", "application/gzip"],
  [".tgz", "application/gzip"],
  [".rds", "application/octet-stream"],
  [".json", "application/json"],
  [".html", "text/html; charset=utf-8"],
  [".txt", "text/plain; charset=utf-8"]
]);

const server = createServer(async (request, response) => {
  const pathname = decodeURIComponent(
    new URL(request.url, "http://127.0.0.1").pathname
  );
  const filePath = resolve(repositoryRoot, `.${pathname}`);
  if (filePath !== repositoryRoot && !filePath.startsWith(`${repositoryRoot}${sep}`)) {
    response.writeHead(403).end();
    return;
  }

  try {
    const contents = await readFile(filePath);
    response.writeHead(200, {
      "access-control-allow-origin": "*",
      "content-type": mimeTypes.get(extname(filePath)) ?? "application/octet-stream"
    });
    response.end(contents);
  } catch {
    response.writeHead(404).end("Not found");
  }
});

await new Promise((resolveListen, rejectListen) => {
  server.once("error", rejectListen);
  server.listen(0, "127.0.0.1", resolveListen);
});

const address = server.address();
const webR = new WebR();

try {
  await webR.init();
  await webR.installPackages("SPmodels", {
    repos: `http://127.0.0.1:${address.port}/`,
    mount: false
  });
  await webR.evalRVoid(`
    library(SPmodels)

    coords <- cbind(x = c(0, 1, 0, 2), y = c(0, 0, 1, 1))
    local_covariance <- array(0, dim = c(nrow(coords), 2, 2))
    for (i in seq_len(nrow(coords))) {
      local_covariance[i, , ] <- diag(1, 2)
    }
    correlation <- sp_nonstationary_cor(
      coords, local_covariance, smoothness = 1.5
    )
    stopifnot(identical(dim(correlation), c(4L, 4L)))
    stopifnot(max(abs(correlation - t(correlation))) < 1e-12)

    dat <- data.frame(
      response = c(1.2, 0.7, 2.1, 1.8, 3.0, 2.6),
      exposure = c(-1, 0, 0.5, 1, 1.5, 2),
      east = c(0, 1, 0, 2, 1, 3),
      north = c(0, 0, 1, 1, 2, 2)
    )
    control <- sp_control(
      smoothness = 0.5,
      spatial_variance = 1.4,
      scale = 1.7,
      nugget_variance = 0.3,
      estimate = FALSE,
      distance_units = "km"
    )
    fit <- sp_fit(
      response ~ exposure,
      data = dat,
      coords = ~ east + north,
      control = control
    )
    prediction <- predict(
      fit,
      newdata = data.frame(exposure = 0.25, east = 0.5, north = 0.25),
      se.fit = TRUE
    )
    stopifnot(inherits(fit, "spmodels_fit"))
    stopifnot(is.finite(as.numeric(logLik(fit))))
    stopifnot(length(prediction$fit) == 1L, is.finite(prediction$se.fit))
    stopifnot(!getFromNamespace(".sp_has_fork_workers", "SPmodels")())
  `);
  console.log(`WebR ${webR.version} / R ${webR.versionR}: SPmodels install and smoke example passed.`);
} finally {
  webR.close();
  await new Promise((resolveClose) => server.close(resolveClose));
}
