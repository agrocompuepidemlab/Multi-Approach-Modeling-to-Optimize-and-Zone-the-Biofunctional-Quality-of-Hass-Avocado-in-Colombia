// Builds a journal-formatted .docx from a simple markup file.
// Usage: node build_docx.js <markup.txt> <out.docx> <figdir> <tabledir> [--supp]
const fs = require('fs');
const path = require('path');
const {
  Document, Packer, Paragraph, TextRun, ImageRun, Table, TableRow, TableCell, WidthType, AlignmentType,
  BorderStyle, ShadingType, Footer, PageNumber, LineNumberRestartFormat, HeadingLevel, PageBreak, LineRuleType,
} = require('docx');

const [, , inFile, outFile, figDir, tabDir] = process.argv;
const SUPP = process.argv.includes('--supp');
const FONT = 'Times New Roman';
const LINE = process.argv.includes('--supp') ? 276 : 480; // double spacing for manuscript

function sizeOfPng(buf) { return { w: buf.readUInt32BE(16), h: buf.readUInt32BE(20) }; }

// inline markup: **bold**, *italic*, ^sup^, ~sub~, \* literal
function runs(text, base = {}) {
  const out = []; let b = false, i = false, sup = false, sub = false, cur = '';
  const flush = () => { if (cur) { out.push(new TextRun({ text: cur, bold: b || base.bold, italics: i !== !!base.italics, superScript: sup, subScript: sub, font: FONT, size: base.size || 24 })); cur = ''; } };
  for (let k = 0; k < text.length; k++) {
    const c = text[k];
    if (c === '\\' && k + 1 < text.length) { cur += text[k + 1]; k++; continue; }
    if (c === '*' && text[k + 1] === '*') { flush(); b = !b; k++; continue; }
    if (c === '*') { flush(); i = !i; continue; }
    if (c === '^') { flush(); sup = !sup; continue; }
    if (c === '~') { flush(); sub = !sub; continue; }
    cur += c;
  }
  flush(); return out;
}
const para = (text, opt = {}) => new Paragraph({
  children: runs(text, opt.run || {}), alignment: opt.align || AlignmentType.JUSTIFIED,
  spacing: { line: opt.line || LINE, lineRule: LineRuleType.AUTO, after: opt.after ?? 0, before: opt.before ?? 0 },
  indent: opt.indent, heading: opt.heading, keepNext: opt.keepNext,
});

function parseCSV(txt) {
  const rows = []; let row = [], f = '', q = false;
  for (let k = 0; k < txt.length; k++) {
    const c = txt[k];
    if (q) { if (c === '"' && txt[k + 1] === '"') { f += '"'; k++; } else if (c === '"') q = false; else f += c; }
    else if (c === '"') q = true; else if (c === ',') { row.push(f); f = ''; }
    else if (c === '\n') { row.push(f); rows.push(row); row = []; f = ''; } else if (c !== '\r') f += c;
  }
  if (f || row.length) { row.push(f); rows.push(row); }
  return rows;
}
const HEADERS = {
  'Table1.csv': ['Trait', 'North', 'East', 'Southwest', 'Main crop', 'Off-season', '*P* subregion', '*P* harvest'],
  'Table2.csv': ['Trait', 'Deviance explained (%)', 'Harvest coef. (*P*)', 's(*T*) edf (*P*)', 's(pH) edf (*P*)', 's(DM) edf (*P*)', '*Q*^2^ LOFO GAM / MLR', 'y-rand. *P* GAM / MLR'],
};
const WIDTHS = { 'Table1_best_models.csv': [1500,1400,2700,800,800,650,750,760], 'TableS1_variables.csv': [1300,3800,2200,2060], 'TableS2_all_models.csv': [2100,700,900,1500,1050,1050,1000,1060], 'TableS3_corrections.csv': [700,2000,1500,5160], 'TableS4_zoning_area.csv': [2560,1700,1700,1700,1700], 'TableS1.csv': [1100,2900,1700,1500,900,1260], 'TableS4.csv': [1400,2000,2000,2000,980,980], 'TableS7.csv': [1500,3300,2200,2360], 'TableS8.csv': [1300,900,2000,4060,1100], 'TableS3.csv': [1500,1700,1100,1300,1300,1160,1300], 'Table1.csv': [1600,1248,1248,1248,1248,1248,760,760], 'Table2.csv': [1500,800,1250,1150,1150,1150,1180,1180] };
function table(csvName, widthsOpt) {
  widthsOpt = widthsOpt || WIDTHS[csvName];
  const rows = parseCSV(fs.readFileSync(path.join(tabDir, csvName), 'utf8').trim());
  const header = HEADERS[csvName] || rows[0];
  const body = rows.slice(1);
  const n = header.length; const total = 9360;
  const widths = widthsOpt || (() => { const first = Math.round(total * 0.22); const rest = Math.floor((total - first) / (n - 1)); const w = [first, ...Array(n - 1).fill(rest)]; w[n - 1] += total - w.reduce((a, b) => a + b, 0); return w; })();
  const border = { style: BorderStyle.SINGLE, size: 4, color: '000000' };
  const none = { style: BorderStyle.NONE, size: 0, color: 'FFFFFF' };
  const mk = (txt, w, isHead, isLast) => new TableCell({
    width: { size: w, type: WidthType.DXA },
    borders: { top: isHead ? border : none, bottom: (isHead || isLast) ? border : none, left: none, right: none },
    shading: isHead ? { type: ShadingType.CLEAR, color: 'auto', fill: 'E7E6E6' } : undefined,
    margins: { top: 30, bottom: 30, left: 30, right: 30 },
    children: [new Paragraph({ children: runs(txt, { size: 16, bold: isHead }), alignment: AlignmentType.LEFT, spacing: { line: 240, lineRule: LineRuleType.AUTO } })],
  });
  const trs = [new TableRow({ tableHeader: true, children: header.map((h, j) => mk(h, widths[j], true, false)) })];
  body.forEach((r, ri) => trs.push(new TableRow({ children: r.map((c, j) => mk(c, widths[j], false, ri === body.length - 1)) })));
  return new Table({ width: { size: total, type: WidthType.DXA }, columnWidths: widths, rows: trs });
}
function figure(file) {
  const buf = fs.readFileSync(path.join(figDir, file)); const { w, h } = sizeOfPng(buf);
  const maxW = 624, maxH = 760; let W = maxW, H = Math.round(h * maxW / w); if (H > maxH) { H = maxH; W = Math.round(w * maxH / h); }
  return new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 120, after: 60, line: 240, lineRule: LineRuleType.AUTO }, keepNext: true,
    children: [new ImageRun({ type: 'png', data: buf, transformation: { width: W, height: H } })] });
}

const children = [];
const lines = fs.readFileSync(inFile, 'utf8').split('\n');
for (const raw of lines) {
  if (!raw.trim()) continue;
  const m = raw.match(/^([A-Z0-9]+):\s?(.*)$/); if (!m) { children.push(para(raw)); continue; }
  const [, tag, t] = m;
  switch (tag) {
    case 'TITLE': children.push(para(t, { align: AlignmentType.CENTER, run: { bold: true, size: 28 }, after: 120 })); break;
    case 'STITLE': children.push(para(t, { align: AlignmentType.CENTER, run: { italics: true, size: 24 }, after: 240 })); break;
    case 'AUTH': children.push(para(t, { align: AlignmentType.JUSTIFIED, after: 120 })); break;
    case 'AFF': children.push(para(t, { align: AlignmentType.LEFT, run: { size: 20 } })); break;
    case 'NOTE': children.push(para(t, { align: AlignmentType.LEFT, run: { size: 20, italics: true }, after: 120 })); break;
    case 'H1': children.push(para(t, { align: AlignmentType.LEFT, run: { bold: true }, before: 240, keepNext: true, heading: HeadingLevel.HEADING_1 })); break;
    case 'H2': children.push(para(t, { align: AlignmentType.LEFT, run: { bold: true, italics: true }, before: 120, keepNext: true, heading: HeadingLevel.HEADING_2 })); break;
    case 'P': children.push(para(t, { indent: { firstLine: 0 }, after: 0 })); break;
    case 'KW': children.push(para(t, { after: 120 })); break;
    case 'EQ': children.push(para(t, { align: AlignmentType.CENTER })); break;
    case 'REF': children.push(para(t, { align: AlignmentType.LEFT, indent: { left: 720, hanging: 720 } })); break;
    case 'PB': children.push(new Paragraph({ children: [new PageBreak()] })); break;
    case 'FIG': { const [file, cap] = t.split('|'); children.push(figure(file)); children.push(para(cap, { align: AlignmentType.JUSTIFIED, after: 240 })); break; }
    case 'TAB': { const [file, cap, foot] = t.split('|'); children.push(para(cap, { align: AlignmentType.JUSTIFIED, keepNext: true, after: 60 }));
      children.push(table(file)); if (foot) children.push(para(foot, { run: { size: 20 }, line: 360, after: 240, before: 60 })); break; }
    default: children.push(para(t));
  }
}

const doc = new Document({
  styles: { default: { document: { run: { font: FONT, size: 24 } } },
    paragraphStyles: [
      { id: 'Heading1', name: 'Heading 1', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: { font: FONT, size: 24, bold: true, color: '000000' }, paragraph: { outlineLevel: 0 } },
      { id: 'Heading2', name: 'Heading 2', basedOn: 'Normal', next: 'Normal', quickFormat: true, run: { font: FONT, size: 24, bold: true, italics: true, color: '000000' }, paragraph: { outlineLevel: 1 } },
    ] },
  sections: [{
    properties: {
      page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } },
      lineNumbers: SUPP ? undefined : { countBy: 1, restart: LineNumberRestartFormat.CONTINUOUS, distance: 360 },
    },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 20 })] })] }) },
    children,
  }],
});
Packer.toBuffer(doc).then(b => { fs.writeFileSync(outFile, b); console.log('wrote', outFile); });
