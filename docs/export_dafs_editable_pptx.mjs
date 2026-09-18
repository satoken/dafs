#!/usr/bin/env node

/**
 * Export the project-local HTML deck as a native, editable PPTX.
 *
 * The exporter is intentionally kept as a small wrapper around dom-to-pptx:
 * it measures the already-rendered DOM and maps text, boxes, tables and SVG
 * elements to PowerPoint objects instead of taking a screenshot of a slide.
 *
 * The local Firefox installation is a Snap build whose launcher is not
 * usable by Puppeteer in this environment.  We therefore patch only the
 * package's browser-discovery list in a temporary copy and keep the package
 * itself untouched.
 */

import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import { pathToFileURL } from 'node:url';

const repoRoot = path.resolve(path.dirname(new URL(import.meta.url).pathname), '..');
const htmlPath = path.join(repoRoot, 'docs', 'dafs_linearization_slides.html');
const outputPath = path.join(repoRoot, 'docs', 'dafs_linearization_slides_editable.pptx');
const packageRoot = process.env.DOM_TO_PPTX_ROOT || '/tmp/dafs-dom-to-pptx/node_modules/dom-to-pptx';
const nodeExporterPath = path.join(packageRoot, 'dist', 'dom-to-pptx-node.mjs');
const directFirefox = '/snap/firefox/current/usr/lib/firefox/firefox';

if (!fs.existsSync(htmlPath)) throw new Error(`HTML deck not found: ${htmlPath}`);
if (!fs.existsSync(nodeExporterPath)) {
  throw new Error(
    `dom-to-pptx is not installed at ${packageRoot}. Set DOM_TO_PPTX_ROOT or install dom-to-pptx first.`
  );
}

const originalExporter = fs.readFileSync(nodeExporterPath, 'utf8');
const patchedExporter = originalExporter
  // Puppeteer 25 tries to resolve its bundled Chrome revision before the
  // package's system-browser fallback.  This installation intentionally has
  // no Puppeteer-managed browser, so skip that probe.
  .replace('const p = puppeteer.executablePath();', 'const p = null;')
  .replace(
    "paths: ['/usr/bin/firefox', '/snap/bin/firefox', '/usr/lib/firefox/firefox'],",
    `paths: ['${directFirefox}', '/usr/lib/firefox/firefox'],`
  );
if (patchedExporter === originalExporter) {
  throw new Error('Could not patch dom-to-pptx Firefox discovery path; package layout may have changed.');
}

const tempExporterPath = path.join(packageRoot, '.dafs-dom-to-pptx-node.mjs');
fs.writeFileSync(tempExporterPath, patchedExporter);

try {
  const { exportHtmlToPptx } = await import(`${pathToFileURL(tempExporterPath).href}?dafs=${Date.now()}`);
  const source = fs.readFileSync(htmlPath, 'utf8');

  // The interactive HTML intentionally shows one slide at a time.  The
  // exporter receives an inline copy with all slide roots visible so every
  // .slide becomes a separate PPTX slide.
  const exportMode = `
<script data-dafs-export="true">
  document.querySelectorAll('.slide').forEach((slide) => {
    slide.style.display = 'block';
    slide.style.margin = '0 auto';
  });
</script>`;
  const exportSource = source.replace('</body>', `${exportMode}\n</body>`);

  const buffer = await exportHtmlToPptx(exportSource, {
    selector: '.slide',
    browserWidth: 2560,
    browserHeight: 1440,
    pptxOptions: {
      title: 'DAFS 線形化：実装差分と実験結果',
      author: 'DAFS project',
      width: 13.333,
      height: 7.5,
      includePseudoElements: true,
      svgAsVector: true,
    },
  });

  fs.writeFileSync(outputPath, buffer);
  console.log(`Wrote editable PPTX: ${outputPath}`);
} finally {
  fs.rmSync(tempExporterPath, { force: true });
}
