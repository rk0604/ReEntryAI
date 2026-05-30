import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import { existsSync, createReadStream, statSync } from 'node:fs'
import { resolve, join } from 'node:path'

// Where the simulation writes its outputs. The dashboard serves these
// files LIVE — re-run the sim and refresh the browser, no copy needed.
const REVISION_DIR = resolve(__dirname, '../StochasticEntrySim/revision_v1')

const MIME = {
  json: 'application/json; charset=utf-8',
  csv:  'text/csv; charset=utf-8',
  png:  'image/png',
  jpg:  'image/jpeg',
  jpeg: 'image/jpeg',
  svg:  'image/svg+xml',
}

/**
 * Middleware: serve /data/<file> directly from revision_v1/<file>.
 * No caching, so changes are picked up on the very next request.
 */
function serveLiveRevisionData() {
  return {
    name: 'serve-revision-v1-live',
    configureServer(server) {
      server.middlewares.use('/data', (req, res, next) => {
        const reqPath = (req.url || '/').split('?')[0]
        const filePath = join(REVISION_DIR, reqPath)
        if (!filePath.startsWith(REVISION_DIR)) {
          return next()  // prevent ../ escape
        }
        if (existsSync(filePath) && statSync(filePath).isFile()) {
          const ext = filePath.split('.').pop().toLowerCase()
          res.setHeader('Content-Type', MIME[ext] || 'application/octet-stream')
          res.setHeader('Cache-Control', 'no-store')
          createReadStream(filePath).pipe(res)
        } else {
          next()
        }
      })
    },
  }
}

export default defineConfig({
  plugins: [
    react(),
    serveLiveRevisionData(),
  ],
})
